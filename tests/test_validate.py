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
