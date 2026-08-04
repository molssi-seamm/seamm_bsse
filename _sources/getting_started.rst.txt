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
