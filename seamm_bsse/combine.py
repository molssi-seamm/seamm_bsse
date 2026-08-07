"""Combine the results of the 2N + 1 counterpoise jobs into a corrected total
energy and gradient.

Padding-only, per the architecture doc: the ``cluster`` and
``fragment-in-cluster`` jobs are already full-cluster-length (ghosts keep
their place in the atom list), so only the ``fragment-alone`` jobs -- genuine
subsets, by construction -- need to be padded into the full-cluster-length
gradient before summing.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

from .job_specs import CLUSTER, FRAGMENT_ALONE, JobSpec

#: Net-force tolerance (E_h/bohr) above which a CP-corrected gradient is
#: treated as corrupted -- a translational-invariance guard. A converged
#: CP-corrected gradient of an isolated cluster always sums to ~0 net force;
#: ORCA (RIJCOSX) can silently corrupt the far-ghost Pulay force on a
#: ``fragment-in-cluster`` job while the energy stays fine (no warning, no
#: linear dependence flagged), and this sum is the only reliable detector.
#:
#: Calibrated (2026-08-07) against a real 40-point Na+/Cl- R-scan (job 3867):
#: across the 37 healthy points the combined net force never exceeds
#: 9.6e-5 E_h/bohr; the 3 anomalous points (all traced to exactly this
#: RIJCOSX/ghost-centre noise in a ``fragment-in-cluster`` job, confirmed by
#: hand-recombining the raw per-job gradients -- NOT a wrong SCF branch,
#: which was ruled out by checking Mulliken charges/orbital degeneracy on
#: the same points) sit at 9.5e-4 to 1.3e-3 -- a clean >10x gap. This value
#: is the geometric mean of that gap. The prior value (0.02) was ~200x too
#: loose to catch this class of point at all. See the workspace memory
#: ``nacl-curve-crossing-investigation.md`` for the full investigation
#: (including why a real, distinct, and much smaller symmetry-breaking
#: effect near the ionic/covalent curve crossing does NOT trip this guard --
#: correctly, since the true BSSE correction is already negligible wherever
#: that effect shows up).
#:
#: Where this fires, the physical BSSE correction to the forces is
#: negligible anyway (see ``gradient_correction_magnitude`` on the returned
#: :class:`CPResult` to check this directly for a given point), so falling
#: back to the uncorrected cluster gradient is safe.
DEFAULT_GRADIENT_TOLERANCE = 3e-4


@dataclass(frozen=True)
class JobResult:
    """One job's outcome, as an engine step reports it.

    ``gradient``, if given, must be in the same order as the owning
    ``JobSpec.atom_indices`` (i.e. ``gradient[k]`` is the force on the atom
    at ``spec.atom_indices[k]``), in E_h/bohr; ``None`` for an energy-only
    run.
    """

    energy: float
    gradient: Optional[Sequence[Sequence[float]]] = None


@dataclass(frozen=True)
class CPResult:
    """The counterpoise-corrected outcome for one cluster geometry.

    ``energy`` is the corrected **total** energy of the cluster (same
    absolute scale as the cluster job's raw energy), not an interaction
    energy. ``gradient`` is full-cluster length (E_h/bohr), or ``None`` if
    any job result lacked a gradient. ``bsse_correction`` is
    ``energy - <raw cluster energy>`` (negative: BSSE always raises the
    apparent binding, so the correction lowers the magnitude of binding).
    ``interaction_energy``/``uncorrected_interaction_energy`` are relative to
    the separated fragments (CP-corrected and not, respectively) -- the
    conventional binding-energy reference point, distinct from ``energy``.
    ``gradient_fallback`` is ``True`` when the translational-invariance guard
    fired and ``gradient`` was replaced by the uncorrected cluster gradient
    (``energy`` is unaffected either way). ``net_force`` is the magnitude
    (E_h/bohr) of the net force summed over the combined gradient before any
    fallback -- ``None`` when there is no gradient. ``net_force`` measures
    how *inconsistent* (noisy) the correction is, which is a different
    question from how *large* it is -- ``gradient_correction_magnitude`` is
    the size (E_h/bohr, a full-cluster Frobenius norm) of the BSSE gradient
    correction actually applied (``gradient - <raw cluster gradient>``),
    always computed when a gradient is available (whether or not the
    fallback fires). A caller wondering whether the uncorrected-gradient
    fallback is trustworthy for a given point -- rather than just trusting
    the "large separation" assumption baked into the guard's tolerance --
    can check this directly: a small value means there was little physical
    correction to lose either way; a large value together with
    ``gradient_fallback`` means the fallback is discarding something that
    actually mattered, and the point is worth excluding/rerunning rather
    than silently accepted.
    """

    energy: float
    gradient: Optional[List[List[float]]]
    bsse_correction: float
    interaction_energy: float
    uncorrected_interaction_energy: float
    gradient_fallback: bool = False
    net_force: Optional[float] = None
    gradient_correction_magnitude: Optional[float] = None


def _pad(values, atom_indices, n_atoms):
    """Place a job-local ``[len(atom_indices)][3]`` array into an
    ``[n_atoms][3]`` array of zeros, at the job's global atom indices."""
    if len(values) != len(atom_indices):
        raise ValueError(
            f"Expected {len(atom_indices)} gradient rows (one per atom in "
            f"this job), got {len(values)}."
        )
    padded = [[0.0, 0.0, 0.0] for _ in range(n_atoms)]
    for local_index, global_index in enumerate(atom_indices):
        padded[global_index] = list(values[local_index])
    return padded


def combine(
    specs: Sequence[JobSpec],
    results: Dict[str, JobResult],
    n_atoms: int,
    gradient_tolerance: float = DEFAULT_GRADIENT_TOLERANCE,
) -> CPResult:
    """Assemble the counterpoise-corrected energy/gradient from the 2N + 1
    job results.

    Parameters
    ----------
    specs : Sequence[JobSpec]
        The job specs from :func:`~seamm_bsse.job_specs.generate_job_specs`.
    results : dict[str, JobResult]
        The outcome of each job, keyed by ``JobSpec.label``. Must have an
        entry for every label in ``specs``.
    n_atoms : int
        The number of atoms in the full cluster (the length of the returned
        gradient).
    gradient_tolerance : float
        Net-force tolerance (E_h/bohr) for the translational-invariance
        guard -- see :data:`DEFAULT_GRADIENT_TOLERANCE`.
    """
    by_label = {spec.label: spec for spec in specs}
    missing = [spec.label for spec in specs if spec.label not in results]
    if missing:
        raise KeyError(f"No result given for job(s): {', '.join(missing)}")

    cluster_spec = next(spec for spec in specs if spec.kind == CLUSTER)
    cluster_result = results[cluster_spec.label]
    fragment_labels = [
        spec.fragment_label for spec in specs if spec.kind == FRAGMENT_ALONE
    ]

    have_gradients = all(results[spec.label].gradient is not None for spec in specs)

    e_in_cluster_sum = 0.0
    e_alone_sum = 0.0
    cluster_gradient = None
    correction = None
    if have_gradients:
        cluster_gradient = _pad(
            cluster_result.gradient, cluster_spec.atom_indices, n_atoms
        )
        correction = [[0.0, 0.0, 0.0] for _ in range(n_atoms)]

    for fragment_label in fragment_labels:
        in_cluster_spec = by_label[f"{fragment_label}-in-cluster"]
        alone_spec = by_label[f"{fragment_label}-alone"]
        in_cluster_result = results[in_cluster_spec.label]
        alone_result = results[alone_spec.label]

        e_in_cluster_sum += in_cluster_result.energy
        e_alone_sum += alone_result.energy

        if have_gradients:
            g_in_cluster = _pad(
                in_cluster_result.gradient, in_cluster_spec.atom_indices, n_atoms
            )
            g_alone = _pad(alone_result.gradient, alone_spec.atom_indices, n_atoms)
            correction = [
                [rc + (ga - gi) for rc, gi, ga in zip(row_corr, row_i, row_a)]
                for row_corr, row_i, row_a in zip(correction, g_in_cluster, g_alone)
            ]

    gradient_correction_magnitude = None
    if have_gradients:
        gradient = [
            [gc + corr for gc, corr in zip(row_c, row_corr)]
            for row_c, row_corr in zip(cluster_gradient, correction)
        ]
        gradient_correction_magnitude = (
            sum(c * c for row in correction for c in row) ** 0.5
        )
    else:
        gradient = None

    e_cluster = cluster_result.energy
    delta = e_in_cluster_sum - e_alone_sum  # sum of per-fragment BSSE gaps
    energy = e_cluster - delta

    # Translational-invariance guard: a correct CP-corrected gradient of an
    # isolated cluster always sums to ~0 net force (the correction only
    # redistributes force between fragments' ghost/real centres, it adds
    # none). An engine's per-job gradient can come back silently corrupted
    # (ORCA's RIJCOSX on far diffuse ghosts is the confirmed case -- see
    # DEFAULT_GRADIENT_TOLERANCE's docstring), and this net force is the
    # only detector. When it fires, fall back to the uncorrected cluster
    # gradient for the forces -- safe because the physical BSSE force
    # correction is negligible wherever the guard trips (large separation).
    # The energy is unaffected either way.
    gradient_fallback = False
    net_force = None
    if gradient is not None:
        net = [sum(row[c] for row in gradient) for c in range(3)]
        net_force = sum(c * c for c in net) ** 0.5
        if net_force > gradient_tolerance:
            gradient_fallback = True
            gradient = cluster_gradient

    return CPResult(
        energy=energy,
        gradient=gradient,
        bsse_correction=energy - e_cluster,
        interaction_energy=e_cluster - e_in_cluster_sum,
        uncorrected_interaction_energy=e_cluster - e_alone_sum,
        gradient_fallback=gradient_fallback,
        net_force=net_force,
        gradient_correction_magnitude=gradient_correction_magnitude,
    )
