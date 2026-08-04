=======
History
=======

2026.8.4 -- Initial release: N-fragment counterpoise job-spec generation and combination
    * Added ``Fragment``, ``generate_job_specs``, and ``combine``: given a cluster's
      fragment assignment and per-fragment charge, generates the 2N + 1 job specs a
      full-cluster-basis counterpoise (BSSE) correction needs (the full cluster;
      each fragment in the full cluster basis with the rest ghosted; each fragment
      alone in its own basis) and assembles the 2N + 1 results into a
      counterpoise-corrected total energy and gradient.
    * Added ``validate_fragments`` to check a fragment list for overlapping atoms,
      unsupported open-shell fragments, and (optionally) consistency with a known
      cluster charge/multiplicity.
    * Engine-agnostic by design: no quantum chemistry, no engine-specific input or
      output handling. Turning a job spec into an actual calculation (e.g. ORCA's
      ghost-atom notation) and running it stays in each engine's own SEAMM step --
      ``orca_step``'s ``BSSE`` sub-step is the first consumer.
