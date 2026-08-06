"""Tests for seamm_bsse.combine.combine.

Synthetic energies/gradients, hand-derived expected values -- no QM. The N=2
gradient case uses interleaved (non-contiguous) fragment atom indices
specifically so padding is checked by a real local-to-global index mapping,
not just an identity/leading-block coincidence.
"""

import pytest

from seamm_bsse import (
    DEFAULT_GRADIENT_TOLERANCE,
    Fragment,
    JobResult,
    combine,
    generate_job_specs,
)


def test_n2_energy_only_matches_hand_derivation():
    # Matches the algebra orca_step's bsse.py currently hardcodes for N = 2
    # (energy = e_total - (e_fragA - e_monA) - (e_fragB - e_monB)):
    #   corrected = -10.0 - (-0.1) - (-0.1) = -9.8
    #   bsse_correction = corrected - e_total = 0.2
    #   uncorrected_interaction = e_total - e_monA - e_monB = -0.4
    #   corrected_interaction = e_total - e_fragA - e_fragB = -0.2
    fragments = [
        Fragment(label="A", atom_indices=[0, 1]),
        Fragment(label="B", atom_indices=[2, 3]),
    ]
    specs = generate_job_specs(fragments)
    results = {
        "cluster": JobResult(energy=-10.0),
        "A-in-cluster": JobResult(energy=-4.5),
        "B-in-cluster": JobResult(energy=-5.3),
        "A-alone": JobResult(energy=-4.4),
        "B-alone": JobResult(energy=-5.2),
    }

    result = combine(specs, results, n_atoms=4)

    assert result.energy == pytest.approx(-9.8)
    assert result.bsse_correction == pytest.approx(0.2)
    assert result.uncorrected_interaction_energy == pytest.approx(-0.4)
    assert result.interaction_energy == pytest.approx(-0.2)
    assert result.gradient is None


def test_n2_gradient_padding_interleaved_fragments():
    # A = atoms {0, 2}, B = atoms {1, 3} -- interleaved, so a job-local index
    # does not equal its global index (e.g. A-alone's local index 1 is
    # global atom 2), the real test of _pad's index mapping.
    fragments = [
        Fragment(label="A", atom_indices=[0, 2]),
        Fragment(label="B", atom_indices=[1, 3]),
    ]
    specs = generate_job_specs(fragments)

    def row(x):
        return [x, 0.0, 0.0]

    results = {
        "cluster": JobResult(
            energy=-10.0,
            gradient=[row(1.0), row(2.0), row(3.0), row(4.0)],
        ),
        # atom_indices = (0, 1, 2, 3); real atoms 0, 2; ghosts 1, 3.
        "A-in-cluster": JobResult(
            energy=-4.5,
            gradient=[row(1.1), row(0.02), row(3.1), row(0.02)],
        ),
        # real atoms 1, 3; ghosts 0, 2.
        "B-in-cluster": JobResult(
            energy=-5.3,
            gradient=[row(0.03), row(2.1), row(0.03), row(4.1)],
        ),
        # atom_indices = (0, 2) -- local index 1 is global atom 2.
        "A-alone": JobResult(energy=-4.4, gradient=[row(1.05), row(3.05)]),
        # atom_indices = (1, 3) -- local index 1 is global atom 3.
        "B-alone": JobResult(energy=-5.2, gradient=[row(2.05), row(4.05)]),
    }

    # This synthetic gradient is only built to exercise _pad's index mapping
    # and is not translationally invariant (net force != 0 by construction),
    # so disable the unrelated net-force guard here rather than fight it.
    result = combine(specs, results, n_atoms=4, gradient_tolerance=float("inf"))

    assert result.energy == pytest.approx(-9.8)
    assert result.gradient_fallback is False
    expected_x = [0.92, 1.93, 2.92, 3.93]
    for atom_index, expected in enumerate(expected_x):
        assert result.gradient[atom_index][0] == pytest.approx(expected)
        assert result.gradient[atom_index][1] == pytest.approx(0.0)
        assert result.gradient[atom_index][2] == pytest.approx(0.0)


def test_n3_energy_only_matches_hand_derivation():
    fragments = [
        Fragment(label="A", atom_indices=[0]),
        Fragment(label="B", atom_indices=[1]),
        Fragment(label="C", atom_indices=[2]),
    ]
    specs = generate_job_specs(fragments)
    results = {
        "cluster": JobResult(energy=-30.0),
        "A-in-cluster": JobResult(energy=-10.05),
        "B-in-cluster": JobResult(energy=-10.10),
        "C-in-cluster": JobResult(energy=-10.15),
        "A-alone": JobResult(energy=-10.00),
        "B-alone": JobResult(energy=-10.02),
        "C-alone": JobResult(energy=-10.03),
    }

    result = combine(specs, results, n_atoms=3)

    assert result.energy == pytest.approx(-29.75)
    assert result.bsse_correction == pytest.approx(0.25)
    assert result.uncorrected_interaction_energy == pytest.approx(0.05)
    assert result.interaction_energy == pytest.approx(0.30)


def test_net_force_within_tolerance_no_fallback():
    # A clean (translationally-invariant) N=2 gradient: the combined forces
    # sum to ~0, so the guard must not touch the result.
    fragments = [
        Fragment(label="A", atom_indices=[0]),
        Fragment(label="B", atom_indices=[1]),
    ]
    specs = generate_job_specs(fragments)

    def row(x):
        return [x, 0.0, 0.0]

    results = {
        "cluster": JobResult(energy=-10.0, gradient=[row(1.0), row(-1.0)]),
        "A-in-cluster": JobResult(energy=-4.5, gradient=[row(1.0), row(0.0)]),
        "B-in-cluster": JobResult(energy=-5.3, gradient=[row(0.0), row(-1.0)]),
        "A-alone": JobResult(energy=-4.4, gradient=[row(1.0)]),
        "B-alone": JobResult(energy=-5.2, gradient=[row(-1.0)]),
    }

    result = combine(specs, results, n_atoms=2)

    assert result.gradient_fallback is False
    assert result.net_force == pytest.approx(0.0)
    # in_cluster == alone for both fragments here, so no correction at all.
    assert result.gradient[0][0] == pytest.approx(1.0)
    assert result.gradient[1][0] == pytest.approx(-1.0)


def test_corrupt_ghost_gradient_falls_back_to_cluster_gradient():
    # The ORCA COSX-on-far-ghosts failure mode: ONE fragment's ghost-centre
    # Pulay force blows up (here A-in-cluster's ghost, atom 1) while
    # everything else stays clean -- a one-sided corruption that breaks
    # translational invariance (unlike a symmetric error, which would
    # cancel). The guard must discard it and use the (clean) cluster
    # gradient for the forces, while leaving the CP-corrected energy
    # untouched.
    fragments = [
        Fragment(label="A", atom_indices=[0]),
        Fragment(label="B", atom_indices=[1]),
    ]
    specs = generate_job_specs(fragments)

    def row(x):
        return [x, 0.0, 0.0]

    results = {
        "cluster": JobResult(energy=-10.0, gradient=[row(0.09), row(-0.09)]),
        "A-in-cluster": JobResult(energy=-4.5, gradient=[row(0.1), row(300.0)]),
        "B-in-cluster": JobResult(energy=-5.3, gradient=[row(0.0), row(-0.1)]),
        "A-alone": JobResult(energy=-4.4, gradient=[row(0.1)]),
        "B-alone": JobResult(energy=-5.2, gradient=[row(-0.1)]),
    }

    result = combine(specs, results, n_atoms=2)

    assert result.gradient_fallback is True
    assert result.net_force > DEFAULT_GRADIENT_TOLERANCE
    # Fell back to the clean, uncorrected cluster gradient.
    assert result.gradient[0][0] == pytest.approx(0.09)
    assert result.gradient[1][0] == pytest.approx(-0.09)
    # Energy is still fully CP-corrected -- only the forces fell back.
    assert result.energy == pytest.approx(-9.8)


def test_gradient_tolerance_is_overridable():
    fragments = [
        Fragment(label="A", atom_indices=[0]),
        Fragment(label="B", atom_indices=[1]),
    ]
    specs = generate_job_specs(fragments)

    def row(x):
        return [x, 0.0, 0.0]

    # A small, deliberate net-force residual that is within the default
    # tolerance but not a stricter, caller-supplied one.
    results = {
        "cluster": JobResult(energy=-10.0, gradient=[row(1.0), row(-1.0)]),
        "A-in-cluster": JobResult(energy=-4.5, gradient=[row(1.0), row(0.0)]),
        "B-in-cluster": JobResult(energy=-5.3, gradient=[row(0.0), row(-0.99)]),
        "A-alone": JobResult(energy=-4.4, gradient=[row(1.0)]),
        "B-alone": JobResult(energy=-5.2, gradient=[row(-1.0)]),
    }

    lenient = combine(specs, results, n_atoms=2)
    assert lenient.gradient_fallback is False

    strict = combine(specs, results, n_atoms=2, gradient_tolerance=1e-6)
    assert strict.gradient_fallback is True
    assert strict.gradient[0][0] == pytest.approx(1.0)
    assert strict.gradient[1][0] == pytest.approx(-1.0)


def test_missing_result_raises_keyerror():
    fragments = [
        Fragment(label="A", atom_indices=[0]),
        Fragment(label="B", atom_indices=[1]),
    ]
    specs = generate_job_specs(fragments)
    results = {"cluster": JobResult(energy=-1.0)}

    with pytest.raises(KeyError):
        combine(specs, results, n_atoms=2)


def test_one_missing_gradient_falls_back_to_energy_only():
    fragments = [
        Fragment(label="A", atom_indices=[0]),
        Fragment(label="B", atom_indices=[1]),
    ]
    specs = generate_job_specs(fragments)
    results = {
        "cluster": JobResult(energy=-10.0, gradient=[[0.0, 0.0, 0.0]] * 2),
        "A-in-cluster": JobResult(energy=-4.5, gradient=[[0.0, 0.0, 0.0]] * 2),
        "B-in-cluster": JobResult(energy=-5.3, gradient=[[0.0, 0.0, 0.0]] * 2),
        "A-alone": JobResult(energy=-4.4, gradient=None),  # missing
        "B-alone": JobResult(energy=-5.2, gradient=[[0.0, 0.0, 0.0]]),
    }

    result = combine(specs, results, n_atoms=2)

    assert result.gradient is None
    assert result.energy == pytest.approx(-9.8)
