seamm_bsse
==========

**Status: design only -- not yet implemented.** This package is scaffolded far
enough to hold its design doc; `Fragment`/`JobSpec`/`generate_job_specs`/
`combine` do not exist yet.

An engine-agnostic library for N-fragment, charge-aware counterpoise (BSSE)
corrections: given a cluster + fragment assignment + per-fragment
charge/multiplicity, generate the 2N + 1 job specs (which atoms are real,
ghost, or absent, at what charge/multiplicity, per job) and combine the
2N + 1 results into a CP-corrected total energy and gradient. No quantum
chemistry, no engine-specific input/output handling -- that stays in each
engine's own SEAMM step (`orca_step` first).

See `docs/developer_guide/campaigns/2026-08-03/bsse_architecture.rst` for the
full design: the physics, the confirmed ORCA ghost-indexing detail that makes
`combine()` padding-only, what changes in `orca_step`, and the milestone
sequence.

Motivating plan: `~/Sites/mlff-training/2026-07-28_ion-bsse-plan/`.
