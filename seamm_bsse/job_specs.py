"""Generate the 2N + 1 job specs a full-cluster-basis, N-fragment counterpoise
correction needs, from a list of :class:`~seamm_bsse.fragment.Fragment`.

See the architecture doc
(``docs/developer_guide/campaigns/2026-08-03/bsse_architecture.rst``) for the
physics and the confirmed ORCA ghost-indexing property this shape relies on.
"""

from dataclasses import dataclass
from typing import FrozenSet, List, Optional, Sequence, Tuple

from .fragment import Fragment
from .validate import validate_fragments

#: The three kinds of job a counterpoise correction runs.
CLUSTER = "cluster"
FRAGMENT_IN_CLUSTER = "fragment-in-cluster"
FRAGMENT_ALONE = "fragment-alone"


@dataclass(frozen=True)
class JobSpec:
    """One of the 2N + 1 calculations a counterpoise correction needs.

    Parameters
    ----------
    label : str
        A unique identifier for this job, e.g. ``"cluster"``,
        ``"A-in-cluster"``, ``"A-alone"``.
    kind : str
        One of :data:`CLUSTER`, :data:`FRAGMENT_IN_CLUSTER`,
        :data:`FRAGMENT_ALONE`.
    fragment_label : str | None
        The fragment this job's *real* atoms belong to (the fragment being
        evaluated); ``None`` for the ``"cluster"`` job, which has no single
        owning fragment.
    atom_indices : tuple[int, ...]
        The global (full-cluster) atom indices present in this job, in the
        order a caller's result (energy/gradient) must be returned in. For
        ``CLUSTER`` and ``FRAGMENT_IN_CLUSTER`` jobs this is every atom in
        the cluster, in cluster order (ghosts included) -- see the
        architecture doc's confirmed-for-ORCA indexing property. For
        ``FRAGMENT_ALONE`` jobs it is only that fragment's atoms.
    ghost_indices : frozenset[int]
        The subset of ``atom_indices`` that are ghost centres (basis
        functions only, no nucleus/electrons) in this job.
    charge : int
        The charge to run this job at.
    multiplicity : int
        The spin multiplicity to run this job at.
    """

    label: str
    kind: str
    fragment_label: Optional[str]
    atom_indices: Tuple[int, ...]
    ghost_indices: FrozenSet[int]
    charge: int
    multiplicity: int


def generate_job_specs(fragments: Sequence[Fragment]) -> List[JobSpec]:
    """Return the 2N + 1 :class:`JobSpec` for a cluster of N fragments.

    One ``cluster`` job (every atom real), N ``fragment-in-cluster`` jobs
    (fragment *i* real, every other fragment ghosted, full cluster atom
    list), and N ``fragment-alone`` jobs (only fragment *i*'s atoms, its own
    basis). At N = 2 this is exactly the classic 5-calculation
    Boys-Bernardi scheme.
    """
    validate_fragments(fragments)

    all_atoms = tuple(
        sorted(index for fragment in fragments for index in fragment.atom_indices)
    )
    cluster_charge = sum(fragment.charge for fragment in fragments)

    specs = [
        JobSpec(
            label=CLUSTER,
            kind=CLUSTER,
            fragment_label=None,
            atom_indices=all_atoms,
            ghost_indices=frozenset(),
            charge=cluster_charge,
            multiplicity=1,
        )
    ]

    for fragment in fragments:
        ghosts = frozenset(all_atoms) - frozenset(fragment.atom_indices)
        specs.append(
            JobSpec(
                label=f"{fragment.label}-in-cluster",
                kind=FRAGMENT_IN_CLUSTER,
                fragment_label=fragment.label,
                atom_indices=all_atoms,
                ghost_indices=ghosts,
                charge=fragment.charge,
                multiplicity=fragment.multiplicity,
            )
        )

    for fragment in fragments:
        specs.append(
            JobSpec(
                label=f"{fragment.label}-alone",
                kind=FRAGMENT_ALONE,
                fragment_label=fragment.label,
                atom_indices=fragment.atom_indices,
                ghost_indices=frozenset(),
                charge=fragment.charge,
                multiplicity=fragment.multiplicity,
            )
        )

    return specs
