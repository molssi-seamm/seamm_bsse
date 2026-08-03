"""Tests for seamm_bsse.combine.combine.

Synthetic energies/gradients, hand-derived expected values -- no QM. The N=2
gradient case uses interleaved (non-contiguous) fragment atom indices
specifically so padding is checked by a real local-to-global index mapping,
not just an identity/leading-block coincidence.
"""

import pytest

from seamm_bsse import Fragment, JobResult, combine, generate_job_specs


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

    result = combine(specs, results, n_atoms=4)

    assert result.energy == pytest.approx(-9.8)
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
