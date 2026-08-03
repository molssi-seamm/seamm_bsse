"""The ``Fragment`` input to a counterpoise correction."""

from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class Fragment:
    """One fragment of a cluster for an N-fragment counterpoise correction.

    Parameters
    ----------
    label : str
        A short, unique identifier for the fragment (e.g. ``"A"``, ``"Na+"``,
        ``"water_1"``). Used to name the jobs :func:`generate_job_specs`
        produces for this fragment.
    atom_indices : Sequence[int]
        The 0-based indices of this fragment's atoms in the full cluster's
        atom list.
    charge : int
        The fragment's formal charge.
    multiplicity : int
        The fragment's spin multiplicity. Only ``1`` (closed-shell) is
        currently supported -- see the architecture doc's "out of scope"
        section. Kept as a field now for forward compatibility.
    """

    label: str
    atom_indices: Sequence[int]
    charge: int = 0
    multiplicity: int = 1

    def __post_init__(self):
        indices = tuple(sorted(self.atom_indices))
        if len(indices) == 0:
            raise ValueError(f"Fragment {self.label!r} has no atoms.")
        if len(set(indices)) != len(indices):
            raise ValueError(f"Fragment {self.label!r} has duplicate atom indices.")
        if any(i < 0 for i in indices):
            raise ValueError(f"Fragment {self.label!r} has a negative atom index.")
        # Frozen dataclass -- normalize via object.__setattr__.
        object.__setattr__(self, "atom_indices", indices)
