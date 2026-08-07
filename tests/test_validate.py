"""Tests for seamm_bsse.validate.validate_fragments."""

import pytest

from seamm_bsse import Fragment, validate_fragments


def test_valid_fragments_pass():
    fragments = [
        Fragment(label="A", atom_indices=[0, 1]),
        Fragment(label="B", atom_indices=[2, 3]),
    ]
    validate_fragments(fragments)  # no exception


def test_no_fragments_raises():
    with pytest.raises(ValueError):
        validate_fragments([])


def test_overlapping_atoms_raises():
    fragments = [
        Fragment(label="A", atom_indices=[0, 1]),
        Fragment(label="B", atom_indices=[1, 2]),
    ]
    with pytest.raises(ValueError, match="claimed by both"):
        validate_fragments(fragments)


def test_open_shell_fragment_raises():
    fragments = [
        Fragment(label="A", atom_indices=[0], multiplicity=2),
        Fragment(label="B", atom_indices=[1]),
    ]
    with pytest.raises(ValueError, match="multiplicity"):
        validate_fragments(fragments)


def test_cluster_charge_mismatch_raises():
    fragments = [
        Fragment(label="A", atom_indices=[0], charge=1),
        Fragment(label="B", atom_indices=[1], charge=-1),
    ]
    with pytest.raises(ValueError, match="sum to"):
        validate_fragments(fragments, cluster_charge=1)


def test_cluster_charge_match_passes():
    fragments = [
        Fragment(label="A", atom_indices=[0], charge=1),
        Fragment(label="B", atom_indices=[1], charge=-1),
    ]
    validate_fragments(fragments, cluster_charge=0)  # no exception


def test_cluster_multiplicity_gt_one_raises():
    fragments = [Fragment(label="A", atom_indices=[0])]
    with pytest.raises(ValueError, match="multiplicity"):
        validate_fragments(fragments, cluster_multiplicity=3)


def test_odd_electron_fragment_raises():
    # O, H, H (water, 10 electrons neutral) + Na (11 electrons neutral), with
    # the +1 charge mistakenly assigned to the water fragment instead of Na:
    # the cluster charge (+1) still checks out, but water-at-+1 is a 9
    # electron fragment, which cannot be closed-shell.
    atomic_numbers = [8, 1, 1, 11]
    fragments = [
        Fragment(label="1", atom_indices=[0, 1, 2], charge=1),
        Fragment(label="2", atom_indices=[3], charge=0),
    ]
    with pytest.raises(ValueError, match="9 electrons"):
        validate_fragments(fragments, cluster_charge=1, atomic_numbers=atomic_numbers)


def test_even_electron_fragments_pass():
    # Same cluster, charge assigned to the right fragment (Na+).
    atomic_numbers = [8, 1, 1, 11]
    fragments = [
        Fragment(label="1", atom_indices=[0, 1, 2], charge=0),
        Fragment(label="2", atom_indices=[3], charge=1),
    ]
    validate_fragments(
        fragments, cluster_charge=1, atomic_numbers=atomic_numbers
    )  # no exception


def test_no_atomic_numbers_skips_electron_check():
    # Without atomic_numbers, the odd-electron check is not run (existing
    # callers that don't pass it keep their prior behavior).
    fragments = [Fragment(label="A", atom_indices=[0], charge=1)]
    validate_fragments(fragments)  # no exception
