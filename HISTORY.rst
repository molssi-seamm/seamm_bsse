=======
History
=======

2026.8.6 -- Bugfix: guard against a corrupted counterpoise gradient at large separation
    * ``combine()`` now checks that the assembled counterpoise-corrected gradient
      sums to ~zero net force, as it must by translational invariance. Some
      quantum chemistry codes (confirmed for ORCA, at large fragment separation)
      can silently return a corrupted ghost-centre force while the energy stays
      fine. When the net force exceeds a tolerance (``gradient_tolerance``,
      default 0.02 E_h/bohr), ``combine()`` falls back to the uncorrected cluster
      gradient for the forces -- the physical BSSE correction is negligible
      wherever this triggers -- and reports it via two new ``CPResult`` fields,
      ``gradient_fallback`` and ``net_force``. The corrected energy is unaffected.

2026.8.4.1 -- Internal: refresh the published README; no code changes
    * The README published with 2026.8.4 still described the package as
      design-only and not yet implemented -- left over from before
      ``Fragment``/``JobSpec``/``generate_job_specs``/``combine`` were even
      written. Updated to describe what is actually implemented and
      validated, including the ``orca_step`` and ``psi4_step`` BSSE
      sub-steps that consume this library. ``Fragment``/``JobSpec``/
      ``generate_job_specs``/``combine``/``validate_fragments`` are
      unchanged.

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
