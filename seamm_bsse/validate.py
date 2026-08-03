"""Sanity checks on a fragment list, shared by :func:`generate_job_specs` and
callers that want to cross-check fragments against a known cluster
charge/multiplicity (e.g. an engine step checking its fragment table against
the actual system's ``configuration.charge``)."""

from typing import Optional, Sequence

from .fragment import Fragment


def validate_fragments(
    fragments: Sequence[Fragment],
    cluster_charge: Optional[int] = None,
    cluster_multiplicity: Optional[int] = None,
) -> None:
    """Raise ``ValueError`` if ``fragments`` is not a valid partition, or (when
    given) is inconsistent with ``cluster_charge``/``cluster_multiplicity``.

    Checks:

    * at least one fragment;
    * no atom index is claimed by more than one fragment;
    * every fragment is closed-shell (``multiplicity == 1``) -- open-shell
      fragments are not yet supported (see the architecture doc);
    * if given, ``cluster_charge`` equals the sum of the fragment charges;
    * if given, ``cluster_multiplicity`` is ``1`` (the only multiplicity a set
      of closed-shell fragments can combine to, in this phase).
    """
    if len(fragments) == 0:
        raise ValueError("At least one fragment is required.")

    seen = {}
    for fragment in fragments:
        for index in fragment.atom_indices:
            if index in seen:
                raise ValueError(
                    f"Atom {index} is claimed by both fragment "
                    f"{seen[index]!r} and fragment {fragment.label!r}."
                )
            seen[index] = fragment.label

    for fragment in fragments:
        if fragment.multiplicity != 1:
            raise ValueError(
                f"Fragment {fragment.label!r} has multiplicity "
                f"{fragment.multiplicity}; only closed-shell (multiplicity 1) "
                "fragments are supported."
            )

    if cluster_charge is not None:
        total_charge = sum(fragment.charge for fragment in fragments)
        if total_charge != cluster_charge:
            raise ValueError(
                f"The fragment charges sum to {total_charge}, but the cluster "
                f"charge is {cluster_charge}."
            )

    if cluster_multiplicity is not None and cluster_multiplicity != 1:
        raise ValueError(
            f"Cluster multiplicity {cluster_multiplicity} is not supported; "
            "only closed-shell clusters (multiplicity 1) are, in this phase."
        )
