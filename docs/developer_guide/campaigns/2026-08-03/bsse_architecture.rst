An N-fragment, charge-aware counterpoise architecture
=======================================================

:Status: Draft -- design agreed, implementation not started
:Companion: ``orca_step``'s Phase-1 ``BSSE`` sub-step (two fragments, neutral
            singlet only), scoped in
            ``orca_step/docs/developer_guide/campaigns/2026-07-09/bsse_scope.rst``
:Supersedes (in part): ``~/Work/SEAMM/BSSE_general_step_design.rst`` (2026-07-09)
            -- the physics and validation plan there still hold; the
            *architecture* section is superseded by this document. See
            "Relationship to earlier documents" below.
:Motivating plan: ``~/Sites/mlff-training/2026-07-28_ion-bsse-plan/`` -- the
            Na\ :sup:`+`\ /Cl\ :sup:`-`\ /H\ :sub:`2`\ O pilot that requires
            per-fragment charge and N > 2 fragments, which is what forces this
            generalization.

Goal
----

Generalize the existing two-fragment, neutral-singlet ORCA counterpoise (BSSE)
correction to an arbitrary number of fragments with independent per-fragment
charge, **without** building one monolithic, code-reaching-into-everything
driver step (the QCEngine/QCManyBody failure mode, evaluated and rejected in
the ion-BSSE plan). Scope for this phase: **ORCA only.** Psi4 is the planned
second engine and a cross-validation oracle, but is out of scope for the work
described here.

The physics
-----------

Full-cluster-basis site-site counterpoise, the natural generalization of
Boys-Bernardi to N fragments:

.. math::

   E_\text{CP} = E_\text{clust}(\text{clust})
                 - \sum_{i=1}^{N} \big[\, E_i(\text{clust}) - E_i(i) \,\big]

with the analogous gradient. This needs **2N + 1** calculations per geometry:
the full cluster (all atoms real), N calculations of "fragment *i* in the full
cluster basis" (fragment *i* real, every other fragment ghosted), and N
calculations of "fragment *i* alone" (only fragment *i*'s atoms, its own
basis). At N = 2 this is exactly the existing 5-calculation Boys-Bernardi
scheme ``BSSEGradient.cmp`` already runs, so N = 2 is the regression target
(see "Validation", below).

(The motivating lab-notebook page describes this as "N + 2 job specs" in
passing; the correct count, and what this document and the library implement,
is 2N + 1. They coincide only at N = 2.)

Confirmed design detail: ORCA preserves full-cluster indexing
---------------------------------------------------------------

The ion-BSSE plan flagged one thing to confirm before finalizing the
library's job-spec shape: *do the engines' ghost mechanisms keep the full atom
list/order and just flag atoms real/ghost, or do they remove and renumber?*
If the former, ``combine()`` never needs index remapping -- only padding.

**For ORCA this is already confirmed**, by the existing Phase-1 implementation
(``orca_step/orca_step/bsse.py::_ghost_xyz``, lines 149-159): it writes *every*
atom of the complex for every job that uses ghosts, flagging ghost centres
with a trailing ``:`` on the element symbol -- no atom is ever removed or
renumbered. Consequently:

* The **cluster** job and the N **fragment-in-cluster-basis** jobs are already
  full-cluster-length, in cluster order. Their energies/gradients drop
  straight into the sum/difference with zero remapping.
* Only the N **fragment-alone** jobs are a genuine subset (by construction --
  that is what "alone, own basis" means): a job of ``3|F_i|`` gradient
  components that must be **padded** into the full-cluster-length vector at
  fragment *i*'s global atom indices, zero elsewhere.

So ``combine()`` is exactly "padding + sum/subtract," as hoped, for ORCA. The
Psi4 side of this question (``Gh(...)`` fragment syntax) is explicitly
deferred, not assumed -- see "Out of scope for this phase."

Architecture
------------

Split into a small engine-agnostic library and a thin, generalized ORCA
sub-step, per the ion-BSSE plan's decision 3.

``seamm_bsse`` (new package)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pure bookkeeping: no QM, no ORCA/Psi4/molsystem imports, unit-testable with
synthetic numbers alone. Structured like ``seamm_thermochemistry`` (plain
``pyproject.toml`` + ``versioningit``, no conda dependency, no entry points,
no Tk) rather than the plug-in cookiecutter template, since it registers
nothing with the flowchart.

* ``fragment.py`` -- ``Fragment(label, atom_indices, charge, multiplicity=1)``.
  A cluster is just a list of ``Fragment``.
* ``job_specs.py`` -- ``generate_job_specs(fragments) -> list[JobSpec]``,
  producing the 2N + 1 specs described above. Each ``JobSpec`` carries: a
  label/kind (``"cluster"`` / ``"fragment-in-cluster"`` / ``"fragment-alone"``),
  the global atom indices present in that job, which of those are ghost,
  and the job's charge/multiplicity.
* ``combine.py`` -- ``combine(specs, results, n_atoms) -> CPResult`` with
  ``energy``, ``gradient`` (full-cluster length), ``bsse_correction``,
  ``interaction_energy``, ``uncorrected_interaction_energy``. Generalizes the
  N = 2 algebra currently hardcoded in ``orca_step/orca_step/bsse.py`` (lines
  294-295 for the energy, 307-308 for the interaction energies) to arbitrary
  N, doing exactly the padding described above.
* ``validate.py`` (small) -- fragment-charge/multiplicity sanity checks (sum of
  fragment charges equals the cluster's total charge; see "Out of scope" for
  what is and is not validated in this phase).

``orca_step`` changes
~~~~~~~~~~~~~~~~~~~~~~

* Extract a low-level ``run_orca_job(keyword_line, xyz_text, charge,
  multiplicity, extra_blocks=None, extra_files=None, make_wfx=False)`` out of
  the existing ``run_orca`` (``orca_base.py``, line 176). ``run_orca`` today
  builds its geometry block from ``self.get_system_configuration()`` and reads
  ``configuration.charge``/``configuration.spin_multiplicity`` directly; a
  BSSE sub-job is a different atom subset at a different charge on every call,
  so the new primitive takes geometry text and charge/multiplicity as
  arguments instead. Reuses the resource/MPI-environment/command-construction
  internals ``run_orca`` and ``run_orca_compound`` already share.
* Generalize ``bsse.py``: fragment definition becomes a list of N fragments
  (auto-detected molecules, or specified atom sets), each with an independent
  charge (multiplicity fixed at 1 in this phase -- see "Out of scope"). GUI:
  a fragment table (default one row per auto-detected molecule, charge 0),
  following the existing "prevent invalid combinations in the GUI" rule for
  SEAMM steps.
* ``_ghost_xyz`` generalizes from "always the full atom list, one ghost
  fragment" to "the atoms named in a ``JobSpec``, with its ghost subset
  flagged" -- the same trailing-``:`` mechanism, just parameterized.
* ``run()`` calls ``seamm_bsse.generate_job_specs``, runs each spec through
  ``run_orca_job`` in its own sub-directory, and calls ``seamm_bsse.combine``
  on the collected results.

``molsystem`` -- explicitly unchanged
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The earlier general-step design (``BSSE_general_step_design.rst``) treated a
first-class ``is_ghost`` column on ``molsystem`` atoms as a "critical
dependency" -- but that was written for a different architecture (one driver
node abstracting engines behind a Model-Chemistry/MDI ``evaluate()``
interface). Phase 1 never needed it: it hand-writes ghost-flagged xyz text
directly in the engine step, sidestepping molsystem entirely. This
architecture keeps doing that. No molsystem schema change is part of this
work; revisit only if/when a future engine's native ghost mechanism genuinely
needs molsystem-level representation.

The Compound-script path: kept as the N = 2 oracle, then retired
--------------------------------------------------------------------

``bssegradient.cmp`` / ``bssenergy.cmp`` (the ORCA *Compound* scripts Phase 1
drives) are not touched by this work, and are not deleted yet. They:

* stay in ``orca_step/orca_step/data/`` and continue to be invocable through
  ``run_orca_compound``;
* serve as the **regression oracle** for milestone M2 below -- the new
  ``seamm_bsse``-driven path must reproduce their energies and gradients
  (to the precision the water/FEC/EC campaigns were decided at) on existing
  neutral-singlet dimers before anything N > 2 or charged is trusted;
* are removed from ``bsse.py``'s production run path only *after* that
  regression gate passes -- this is a deliberate two-step "confirm the new
  path is at least as good, then retire the old one" sequence, not a
  unilateral replacement, precisely because the Compound path is
  already-validated, previously-debugged machinery (it has one already-fixed
  energy-formula bug and one already-fixed COSX ghost-gradient blowup guard
  behind it -- see the workspace memory ``bsse-corrected-gradients.md``).

A known trade-off of retiring it: the Compound script runs all five (or
2N + 1) sub-calculations as one ORCA process, which can reuse SCF guesses
across the closely-related sub-jobs; ``seamm_bsse``'s multi-job orchestration
runs each as an independent ORCA invocation. This is expected to cost some
wall-time at N = 2 in exchange for N-fragment and per-fragment-charge
generality; it is not expected to be measured or decided until the M2
regression numbers (energy/gradient parity *and* wall-time) are in hand.

Out of scope for this phase
----------------------------

* **Psi4.** The ghost-indexing confirmation above is ORCA-only; Psi4's
  ``Gh(...)`` mechanism must be checked independently before ``seamm_bsse``
  is assumed to need no Psi4-specific remapping. Deferred to the milestone
  that adds the Psi4 sub-step.
* **Open-shell / multiplicity > 1 fragments.** The ion pilot
  (Na\ :sup:`+`, Cl\ :sup:`-`, H\ :sub:`2`\ O) is closed-shell throughout;
  ``Fragment.multiplicity`` exists in the data model for forward-compatibility
  but this phase only exercises (and only validates) multiplicity = 1 on every
  fragment.
* **Hessians.** Not part of the ion-BSSE plan's roadmap for this phase; a
  CP-corrected Hessian (finite-difference over 2N + 1 already-expensive
  gradient sets, or 2N + 1 analytic Hessians) is a materially larger cost and
  design question, left for later.
* **Periodic / non-molecular-cluster systems.** Per-atom CP is ill-defined for
  a periodic solid; this architecture is molecular-cluster only, as the
  earlier design doc also recommended.
* **VMFC / many-body-decomposed counterpoise.** Full-cluster-basis CP (this
  document) gives one corrected total energy/gradient, which is what MLFF
  training needs; the combinatorial many-body variant is not part of this
  work.

Sequence / milestones
----------------------

Mirrors the ion-BSSE plan's roadmap (``2026-07-28_ion-bsse-plan``), steps 1-3,
ORCA only:

============  ===========================================================
Milestone     Deliverable
============  ===========================================================
M0            This document; ``seamm_bsse`` package scaffolded (structure
              only -- no logic yet).
M1            ``seamm_bsse``: ``Fragment``/``JobSpec``/``generate_job_specs``/
              ``combine``, unit-tested with synthetic energies/gradients (no
              QM). Acceptance: index-mapping/padding tests pass for N = 2 and
              N = 3 synthetic clusters.
M2            ``orca_step``: ``run_orca_job`` extracted; ``bsse.py`` driven by
              ``seamm_bsse`` at N = 2, neutral singlet. **Acceptance:
              reproduces the existing water/FEC/EC Compound-script CP
              energies and gradients** to the precision those campaigns were
              decided at. Compound path retired from the production run path
              once this passes.
M3            Per-fragment charge. Na\ :sup:`+`\ ···H\ :sub:`2`\ O,
              Cl\ :sup:`-`\ ···H\ :sub:`2`\ O, Na\ :sup:`+`\ ···Cl\ :sup:`-`
              validated against literature binding energies (the ion-BSSE
              plan's two-body work).
============  ===========================================================

N = 3 (Na\ :sup:`+`/Cl\ :sup:`-`/H\ :sub:`2`\ O trimers) and the Psi4 sub-step
are the ion-BSSE plan's next roadmap step, deliberately not scoped here.

Validation plan
----------------

#. **Synthetic index-mapping tests** (``seamm_bsse``, no QM): hand-built
   fragment/atom-index cases with made-up per-job energies/gradients, checked
   against hand-computed ``E_CP``/``\nabla E_CP`` for N = 2 and N = 3.
#. **Reproduce Phase 1 exactly** (M2's acceptance gate): same water/FEC/EC
   dimers, same method/basis, through both the Compound-script path and the
   new ``seamm_bsse``-driven path; energies and gradients must agree to the
   precision those campaigns were decided at.
#. **Finite-difference gradient check** for a 3-fragment case once M3 lands,
   since no existing trusted reference covers N = 3.
#. **Literature binding energies** for the charged two-body cases (M3).

Relationship to earlier documents
------------------------------------

* ``orca_step/docs/developer_guide/campaigns/2026-07-09/bsse_scope.rst`` --
  Phase 1. Unchanged by this work except for its eventual role change
  (production path -> regression oracle -> retired from ``bsse.py``, per
  "The Compound-script path" above). Its scope-boundary list (two fragments,
  one charge/mult for all sub-calculations, analytic-gradient methods only) is
  exactly what this document lifts.
* ``~/Work/SEAMM/BSSE_general_step_design.rst`` (2026-07-09) -- the original
  "general, code-agnostic BSSE step" design. Its physics (the 2N + 1 formula),
  results/metadata naming, and validation-plan items are still correct and
  reused above. **Superseded**: its architecture section (one monolithic
  driver node, engines abstracted behind a Model-Chemistry/MDI ``evaluate()``
  interface, ``molsystem`` ghost-atom column as a hard prerequisite) --
  replaced by the library-plus-thin-per-engine-substep split here, per the
  ion-BSSE plan's decision 3 (rejecting the QCEngine-shaped monolith).
* ``~/Sites/mlff-training/2026-07-28_ion-bsse-plan/`` -- the motivating plan
  (the "why"; not duplicated here). This document is the implementation-side
  companion for its "BSSE architecture" section and roadmap steps 1-3.

Open questions
--------------

* Wall-time cost of N x independent ORCA jobs vs. one Compound process, once
  M2's numbers are in hand -- may (or may not) motivate keeping a
  Compound-script fast path for the common N = 2, neutral-singlet case
  alongside the general one, rather than fully retiring it. Deliberately left
  open until the M2 data exists, rather than decided in advance.
* Whether ``seamm_bsse``'s job-level parallelism (the 2N + 1 jobs of one
  geometry are independent) should be exploited now or deferred to a later
  performance pass; not needed for correctness.
