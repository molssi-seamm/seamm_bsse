Getting Started
===============

``seamm_bsse`` is the engine-agnostic library for N-fragment, charge-aware
counterpoise (BSSE) corrections used across SEAMM's quantum-chemistry steps
(``orca_step`` first). Given a cluster's fragment assignment and per-fragment
charge/multiplicity, it generates the 2N + 1 job specs the Boys--Bernardi
counterpoise correction needs (the full cluster; each fragment in the full
cluster basis with the rest ghosted; each fragment alone in its own basis),
and combines the 2N + 1 results into a counterpoise-corrected total energy
and gradient.

It deliberately contains **no quantum chemistry**: no engine-specific input
generation, no ghost-atom notation, no result parsing. Turning a job spec
into an actual ORCA (or, later, Psi4) calculation and running it stays inside
each engine's own SEAMM step. See the
:doc:`architecture design doc <developer_guide/campaigns/2026-08-03/bsse_architecture>`
for the full rationale.

Installing
----------

The library has no runtime dependencies::

    pip install seamm_bsse

Basic usage
-----------

.. code-block:: python

    from seamm_bsse import Fragment, generate_job_specs, combine

    fragments = [
        Fragment(label="A", atom_indices=[0, 1, 2], charge=0),
        Fragment(label="B", atom_indices=[3, 4, 5], charge=0),
    ]
    specs = generate_job_specs(fragments)

    # Run each spec through an engine (not this package's job) to get
    # (energy, gradient) per spec, keyed by spec.label ...
    results = {spec.label: engine.evaluate(spec) for spec in specs}

    result = combine(specs, results, n_atoms=6)
    # result.energy, result.gradient, result.bsse_correction, ...

Choosing per-fragment charges
------------------------------

``Fragment.charge`` is per-fragment, not automatic -- but the convention
followed by every engine step built on this library (``orca_step``, and
``psi4_step``'s independently-implemented equivalent) is:

* If the caller has an explicit charge for each fragment, use it.
* Otherwise, if the source structure carries a **per-atom formal charge**
  (e.g. read from an SDF/MOL file's ``M  CHG`` records), default each
  fragment's charge to the **sum of its own atoms' formal charge**. This
  covers monatomic ions (Na+, Cl-) and polyatomic ions (NH4+, BF4-) alike,
  without the caller needing to specify anything, and is the single most
  effective way to avoid the fragment-ordering mistake described below.
* Otherwise, default to all-neutral.

Adopting this same convention in a new engine integration is strongly
recommended over inventing a different default (for example a
periodic-table/element-based heuristic) -- it generalizes to any ion,
monatomic or not, using information the structure itself already carries,
rather than guessing chemistry.

:func:`validate_fragments` (called by :func:`generate_job_specs`, and
directly by a caller before that if useful) checks two things regardless of
where the charges came from:

* the fragment charges sum to the known cluster charge, if given;
* (when ``atomic_numbers`` is given) each fragment has an **even** number of
  electrons at its assigned charge -- every fragment must be closed-shell.

The second check exists specifically because the first, alone, misses a
charge assigned to the *wrong* fragment when the swap happens to preserve
the cluster's overall charge -- as it always does between two same-parity
fragments, e.g. two monatomic ions. Without it, that mistake only surfaces
as an opaque SCF failure deep inside one engine sub-job. Always pass
``atomic_numbers`` when you have them.

.. tip::

   The fragment *order* -- which is fragment 1, which is fragment 2, ... --
   is whatever order the caller's ``fragments`` list is in (``auto``
   molecule-detection order, or a caller's own atom-group order). It is easy
   to get backwards when supplying explicit charges by hand; a caller-facing
   step should log each fragment's composition and assigned charge before
   running anything, so a user can catch a mismatch before it burns compute.

Diagnosing a corrupted gradient
---------------------------------

:func:`combine` guards the assembled gradient with a translational-
invariance check (a correct counterpoise-corrected gradient of an isolated
cluster always sums to ~0 net force) and falls back to the uncorrected
cluster gradient when that is violated -- see ``CPResult.gradient_fallback``
and ``CPResult.net_force``.

Whether that fallback is actually *safe* for a given point depends on how
much physics it discards, which is a separate question from how noisy the
input gradients were. ``CPResult.gradient_correction_magnitude`` -- the size
of the BSSE gradient correction that was applied (or would have been,
absent any fallback) -- answers that directly:

* ``net_force`` and ``gradient_correction_magnitude`` close in size means
  the "correction" was mostly noise to begin with -- the fallback loses
  little, and is safe.
* ``gradient_correction_magnitude`` much larger than ``net_force`` means
  real physics is being discarded -- the point is worth flagging for
  exclusion or rerun, not silently accepted.

A caller-facing step should report both numbers, not just whether the
fallback fired, so a user can make this judgment -- see ``orca_step``'s BSSE
sub-step user guide for a worked example of the resulting warning message.
