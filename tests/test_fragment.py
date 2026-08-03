"""Tests for seamm_bsse.fragment.Fragment."""

import pytest

from seamm_bsse import Fragment


def test_basic_construction():
    fragment = Fragment(label="A", atom_indices=[2, 0, 1], charge=-1, multiplicity=1)
    assert fragment.label == "A"
    assert fragment.atom_indices == (0, 1, 2)  # sorted
    assert fragment.charge == -1
    assert fragment.multiplicity == 1


def test_defaults():
    fragment = Fragment(label="A", atom_indices=[0])
    assert fragment.charge == 0
    assert fragment.multiplicity == 1


def test_empty_atom_indices_raises():
    with pytest.raises(ValueError):
        Fragment(label="A", atom_indices=[])


def test_duplicate_atom_indices_raises():
    with pytest.raises(ValueError):
        Fragment(label="A", atom_indices=[0, 1, 1])


def test_negative_atom_index_raises():
    with pytest.raises(ValueError):
        Fragment(label="A", atom_indices=[-1, 0])
