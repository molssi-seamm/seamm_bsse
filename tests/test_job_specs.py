"""Tests for seamm_bsse.job_specs.generate_job_specs.

Structural checks only (atom lists, ghost sets, charges) -- the numeric CP
algebra is covered in test_combine.py.
"""

from seamm_bsse import (
    CLUSTER,
    FRAGMENT_ALONE,
    FRAGMENT_IN_CLUSTER,
    Fragment,
    generate_job_specs,
)


def test_n2_neutral_job_count_and_kinds():
    fragments = [
        Fragment(label="A", atom_indices=[0, 1, 2]),
        Fragment(label="B", atom_indices=[3, 4, 5]),
    ]
    specs = generate_job_specs(fragments)

    # 2N + 1 = 5, matching the classic Boys-Bernardi scheme at N = 2.
    assert len(specs) == 5
    kinds = [spec.kind for spec in specs]
    assert kinds.count(CLUSTER) == 1
    assert kinds.count(FRAGMENT_IN_CLUSTER) == 2
    assert kinds.count(FRAGMENT_ALONE) == 2


def test_n2_atom_lists_and_ghosts():
    fragments = [
        Fragment(label="A", atom_indices=[0, 1, 2]),
        Fragment(label="B", atom_indices=[3, 4, 5]),
    ]
    specs = generate_job_specs(fragments)
    by_label = {spec.label: spec for spec in specs}

    # Cluster and fragment-in-cluster jobs are full-cluster-length, in
    # cluster order -- the confirmed-for-ORCA property that makes combine()
    # padding-only.
    full = (0, 1, 2, 3, 4, 5)
    assert by_label["cluster"].atom_indices == full
    assert by_label["cluster"].ghost_indices == frozenset()
    assert by_label["A-in-cluster"].atom_indices == full
    assert by_label["A-in-cluster"].ghost_indices == frozenset({3, 4, 5})
    assert by_label["B-in-cluster"].atom_indices == full
    assert by_label["B-in-cluster"].ghost_indices == frozenset({0, 1, 2})

    # Fragment-alone jobs are genuine subsets -- only this fragment's atoms,
    # no ghosts.
    assert by_label["A-alone"].atom_indices == (0, 1, 2)
    assert by_label["A-alone"].ghost_indices == frozenset()
    assert by_label["B-alone"].atom_indices == (3, 4, 5)
    assert by_label["B-alone"].ghost_indices == frozenset()


def test_n2_charged_fragments_thread_charge_per_job():
    # Na+ / Cl-, the ion-BSSE pilot's simplest charged case.
    fragments = [
        Fragment(label="Na", atom_indices=[0], charge=1),
        Fragment(label="Cl", atom_indices=[1], charge=-1),
    ]
    specs = generate_job_specs(fragments)
    by_label = {spec.label: spec for spec in specs}

    # The cluster (Na+/Cl- pair) is neutral overall even though the
    # fragments are not.
    assert by_label["cluster"].charge == 0
    assert by_label["Na-in-cluster"].charge == 1
    assert by_label["Na-alone"].charge == 1
    assert by_label["Cl-in-cluster"].charge == -1
    assert by_label["Cl-alone"].charge == -1


def test_n3_job_count_and_ghosts_interleaved_fragments():
    # Non-contiguous fragment atom assignment (interleaved), to check the
    # ghost-set/atom-order logic doesn't assume fragments are contiguous
    # blocks.
    fragments = [
        Fragment(label="A", atom_indices=[0, 2], charge=1),
        Fragment(label="B", atom_indices=[1, 3], charge=-1),
        Fragment(label="C", atom_indices=[4]),
    ]
    specs = generate_job_specs(fragments)

    # 2N + 1 = 7 at N = 3 (not "N + 2" = 5, the lab-notebook's informal
    # phrasing that only coincides with 2N + 1 at N = 2).
    assert len(specs) == 7

    by_label = {spec.label: spec for spec in specs}
    full = (0, 1, 2, 3, 4)
    assert by_label["cluster"].atom_indices == full
    assert by_label["A-in-cluster"].ghost_indices == frozenset({1, 3, 4})
    assert by_label["B-in-cluster"].ghost_indices == frozenset({0, 2, 4})
    assert by_label["C-in-cluster"].ghost_indices == frozenset({0, 1, 2, 3})
    assert by_label["A-alone"].atom_indices == (0, 2)
    assert by_label["B-alone"].atom_indices == (1, 3)
    assert by_label["C-alone"].atom_indices == (4,)
