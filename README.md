seamm_bsse
==========
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/molssi-seamm/seamm_bsse/workflows/CI/badge.svg)](https://github.com/molssi-seamm/seamm_bsse/actions?query=workflow%3ACI)

An engine-agnostic library for N-fragment, charge-aware counterpoise (BSSE)
corrections: given a cluster + fragment assignment + per-fragment charge,
`generate_job_specs` produces the 2N + 1 job specs a full-cluster-basis
correction needs (the full cluster; each fragment in the full cluster basis
with the rest ghosted; each fragment alone in its own basis), and `combine`
assembles the 2N + 1 results into a counterpoise-corrected total energy and
gradient. `validate_fragments` checks a fragment list for overlapping atoms,
unsupported open-shell fragments, and consistency with a known cluster
charge/multiplicity.

No quantum chemistry, no engine-specific input/output handling -- that stays
in each engine's own SEAMM step. `orca_step`'s **BSSE** sub-step is the first
consumer: it drives real ORCA jobs through this library's job specs and
`combine()`, generalizing what used to be a two-fragment, neutral-singlet-only
ORCA *Compound*-script correction to an arbitrary number of fragments with
independent per-fragment charge. Validated against real ORCA on water/FEC/EC
dimers (reproduces the old, trusted Compound-script numbers to SCF-noise
level), charged two-body ion/water pairs (Na+/Cl-, Na+/H2O, Cl-/H2O, agreeing
with approximate literature binding energies to a few kcal/mol), and a real
N = 3 Na+/Cl-/H2O trimer. A Psi4 sub-step is planned next, as a thin wrapper
around Psi4's own native `bsse_type='cp'` driver rather than a second consumer
of this library's job-spec/`combine()` path -- Psi4 already does its own
N-fragment CP orchestration internally, so wrapping it gives an independent
cross-check of the ORCA path instead of the same arithmetic run twice.

See `docs/developer_guide/campaigns/2026-08-03/bsse_architecture.rst` for the
full design: the physics, the confirmed ORCA ghost-indexing detail that makes
`combine()` padding-only, what changed in `orca_step`, the validation results,
and the milestone sequence.

Motivating plan: `~/Sites/mlff-training/2026-07-28_ion-bsse-plan/`.
