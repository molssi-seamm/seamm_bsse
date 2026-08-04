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

So ``combine()`` is exactly "padding + sum/subtract," as hoped, for ORCA.
Psi4's side of this question turned out moot rather than merely deferred:
Psi4's native ``bsse_type='cp'`` driver does its own internal ghosting and
combination and hands back an already-full-cluster-indexed result (confirmed
2026-08-04, see "Psi4 sub-step" below), so the Psi4 sub-step never calls
``seamm_bsse.combine()`` at all.

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

* **Psi4.** Not implemented yet. Its architecture is decided (see "Psi4
  sub-step" below) and differs from ORCA's: a thin wrapper around Psi4's
  native ``bsse_type='cp'`` driver, not a second consumer of
  ``generate_job_specs()``/``combine()``.
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
              ``seamm_bsse``, generalized to N fragments/per-fragment charge.
              **Coding done** (2026-08-03): ``run_orca_job`` extracted from
              ``run_orca``; ``geometry_block`` takes ``atom_indices``/
              ``ghost_atoms``; ``bsse.py``/``bsse_parameters.py``/
              ``tk_bsse.py`` generalized (fragment table -> N groups +
              per-fragment charge, validated via ``seamm_bsse.validate_
              fragments``); ``run()`` loops the 2N + 1 job specs through
              ``run_orca_job`` and calls ``seamm_bsse.combine()``. Verified
              with stubbed-ORCA unit tests (job-charge/ghost-set wiring, and
              the assembled energy against hand-derived numbers) -- 135/135
              orca_step tests green, lint clean. **Water-dimer acceptance gate
              PASSED (2026-08-03, real ORCA 6.1.1, HF/def2-SVP,
              ``H2O-H2O.sdf``)**: the new path's corrected energy, BSSE
              correction, and both interaction energies agree with the old
              Compound-script path to ~7e-9 E_h, the gradient to ~7e-8
              E_h/bohr -- SCF-noise level, far inside any precision the
              water/FEC/EC campaigns were decided at. Script:
              ``orca_step/docs/developer_guide/campaigns/2026-07-09/
              validate_bsse_m2.py`` (built a real molsystem Configuration and
              a minimal-but-real execution harness -- actual ``seamm_exec``
              local executor, actual ORCA -- rather than a full flowchart, so
              both the old Compound-script helpers, still present in
              ``bsse.py`` but unused by ``run()``, and the new path could run
              side by side on the identical geometry/method/basis).
              **FEC and EC legs also PASSED (2026-08-03, same script/method,
              ``FEC-FEC dimer opt PM6-ORG.sdf`` and the first record of
              ``EC_dimers_run_1.sdf``)**: energies agree to ~2-3e-9 E_h,
              gradients to ~5.6e-8 E_h/bohr -- the same SCF-noise level as
              water, on 20-atom systems including F (FEC) and a
              production-sampled (NMS-displaced) EC geometry, not just a
              hand-picked equilibrium one. **M2's acceptance gate is
              complete** for N = 2, neutral-singlet ORCA counterpoise; the
              Compound-script path (``bssegradient.cmp``/``bssenergy.cmp``)
              can now be considered fully superseded for production use (kept
              in ``bsse.py`` as a permanent regression reference, not deleted).
M3            Per-fragment charge. Na\ :sup:`+`\ ···H\ :sub:`2`\ O,
              Cl\ :sup:`-`\ ···H\ :sub:`2`\ O, Na\ :sup:`+`\ ···Cl\ :sup:`-`
              validated against literature binding energies (the ion-BSSE
              plan's two-body work). **PASSED (2026-08-04)**, real ORCA,
              B3LYP-D3BJ/def2-TZVP CP interaction energies vs. approximate
              literature references: Na\ :sup:`+`\ ···Cl\ :sup:`-` well
              minimum ≈ -136 kcal/mol at R ≈ 2.4-2.6 Å (lit. ≈ -133 kcal/mol
              at R\ :sub:`e` ≈ 2.36 Å, Born-Haber from atomic D0/IP/EA),
              decaying smoothly to -64 kcal/mol by 7 Å; Na\ :sup:`+`\ ···H\
              :sub:`2`\ O = -26.1 kcal/mol (lit. ≈ -24); Cl\ :sup:`-`\ ···H\
              :sub:`2`\ O = -15.6 kcal/mol (lit. ≈ -13). All within a few
              kcal/mol on **unoptimized, literature-informed geometries** (no
              relaxation at this level of theory) -- exactly the residual
              expected from geometry, not a wiring defect. Confirms
              per-fragment charge (``seamm_bsse.Fragment``/``fragment
              charges``) produces chemically sane numbers on three distinct
              real charged systems (cation-anion, cation-water, anion-water),
              not just correct wiring on stubbed numbers (M2's unit tests).
              Script: ``orca_step/docs/developer_guide/campaigns/2026-07-09/
              validate_bsse_m3.py`` (``"specified"`` fragments -- these
              configurations are hand-built via ``atoms.append`` with no bond
              table, so ``"auto (molecules)"``/``find_molecules()`` would see
              every atom as its own fragment). **Not done**: geometry
              optimization / systematic angular scans (the ion-BSSE plan's
              own next step for these pairs), and the (12)+(3) trimer.
============  ===========================================================

N = 3 on the ORCA path (real ORCA, not just the synthetic N=3 unit tests in
``seamm_bsse`` itself) **PASSED (2026-08-04)**: a hand-built Na\ :sup:`+`\
···Cl\ :sup:`-`\ ···H\ :sub:`2`\ O trimer (contact ion pair, R = 2.44 Å, plus
a water coordinating Na\ :sup:`+`\  from a different direction than Cl\
:sup:`-`), same B3LYP-D3BJ/def2-TZVP level as the M3 two-body checks. Ran the
full 2N + 1 = 7 jobs correctly (confirmed from the printed per-fragment
charges: Na\ :sup:`+`\ /Cl\ :sup:`-`\ /H\ :sub:`2`\ O = +1/&minus;1/0, and
"Running 7 ORCA job(s)"). CP-corrected interaction energy &minus;150.2
kcal/mol: deeper than the isolated Na\ :sup:`+`\ ···Cl\ :sup:`-`\  pair alone
(&minus;136, from M3) but *less* than the naive pairwise sum of the
independently-validated two-body terms (&minus;136 + &minus;26 ≈ &minus;162)
-- the expected cooperative-saturation non-additivity (Na\ :sup:`+`\ , already
partially satisfied by Cl\ :sup:`-`\ , has less capacity left to bind water),
not a bug. Script:
``orca_step/docs/developer_guide/campaigns/2026-07-09/validate_bsse_n3.py``.

Psi4 sub-step: IMPLEMENTED (2026-08-04)
------------------------------------------

The cross-engine check N = 3 was originally paired with -- the Psi4
sub-step -- is done as a first pass: ``psi4_step``'s own ``BSSE`` sub-step
(``molssi-seamm/psi4_step`` PR `#41
<https://github.com/molssi-seamm/psi4_step/pull/41>`_, draft, not yet
merged). Its architecture differs from the ORCA one in an important way.

**Psi4 has a native N-fragment CP driver; ORCA does not.** That native driver
is the whole reason ``seamm_bsse`` exists for ORCA (SEAMM has to do the
2N + 1 job-spec/``combine()`` bookkeeping by hand because ORCA has nothing
built in). Psi4 does not have that gap: ``energy(...)``/``gradient(...)``
with ``bsse_type='cp'`` on a ``--``-separated, per-fragment-charge molecule
block runs the full N-fragment correction internally and returns the
already-combined, full-cluster-indexed result.

Confirmed empirically (2026-08-04, real Psi4 1.10, local ``seamm-psi4`` conda
env -- not the ``psi4.ini`` Docker default, which looks misconfigured,
pointing at a ``seamm-mopac`` container; unrelated, not touched):

* ``gradient('scf', bsse_type='cp', return_wfn=True)`` on the water dimer
  returned a CP-corrected energy *and* a gradient with exactly 6 rows (the
  full 6-atom cluster) in one call -- no remapping needed on the SEAMM side,
  because there is no combine step on the SEAMM side.
* A charged two-fragment case (Na\ :sup:`+`\ /Cl\ :sup:`-`\ , per-fragment
  ``1 1``/``-1 1`` charge/multiplicity headers in the molecule block) gave a
  CP interaction energy of **-136.0 kcal/mol** -- matching the ORCA
  hand-rolled M3 result (-136 kcal/mol) closely, at a different level of
  theory. Free independent cross-validation, before any Psi4 sub-step code
  exists.

**Decision (confirmed with the user 2026-08-04): the Psi4 sub-step is a thin
wrapper around ``bsse_type='cp'``, not a second consumer of
``seamm_bsse.generate_job_specs()``/``combine()``.** This is also the
*better* validation architecture per decision 6 above: two genuinely
independent implementations (ORCA hand-rolled vs. Psi4-native) checking each
other, rather than the same ``combine()`` arithmetic run twice with a
different QM backend underneath. ``seamm_bsse``'s role for Psi4 shrinks to
``Fragment``/``validate_fragments`` (a shared fragment-definition/charge-
validation layer, for GUI consistency with the ORCA sub-step) -- the
job-spec-generation and ``combine()`` pieces go unused for this engine.

One real implementation wrinkle, resolved during implementation:
``psi4_step``'s existing ``Energy`` sub-step builds geometry from the
*global* current configuration (``system_db.system.configuration`` in
``Psi4._convert_structure``), not a per-node ``get_system_configuration()``
like ORCA's sub-steps. ``BSSE`` builds its own fragment-aware,
``--``-separated, per-fragment charge/multiplicity molecule block instead of
reusing ``_convert_structure``. A second wrinkle, also resolved: unlike ORCA,
``psi4_step`` sub-steps don't drive their own execution -- the main ``Psi4``
node concatenates every sub-step's ``get_input()`` text into *one* shared
script and runs it as a single process, so ``BSSE.get_input()`` writes its
JSON result to an **absolute path** (its own sub-step directory), since the
shared process's cwd is the main node's directory, not each sub-step's own.

Validated end-to-end with real Psi4 1.10 on the water dimer (the same
geometry the ORCA M2 regression used): the whole pipeline works, and the
HF/def2-SVP result cross-checks against the independently-validated ORCA
HF/def2-SVP result to ~2e-4 E\ :sub:`h` -- the expected, healthy level of
agreement between two different QC codes at the same nominal level of
theory (not the ~1e-9 machine-precision level the ORCA-vs-its-own-Compound-
script M2 regression showed, which compares the same code against itself).
Script: ``psi4_step/docs/developer_guide/campaigns/2026-08-04/
validate_psi4_bsse.py``.

**Not yet done**: charged fragments and N = 3 validated through this
sub-step specifically (the native Psi4 driver was separately confirmed on
charged fragments directly via hand-written ``psi4`` scripts, not yet
through ``BSSE.get_input()``/``analyze()``); energy-of-formation support;
advanced SCF convergence controls.

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
