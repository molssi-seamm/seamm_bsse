"""Engine-agnostic N-fragment counterpoise (BSSE) job-spec generation and
result combination for SEAMM.

Given a cluster's fragment assignment and per-fragment charge/multiplicity,
generates the 2N + 1 job specs the full-cluster-basis counterpoise correction
needs, and combines the 2N + 1 results into a corrected total energy and
gradient. No quantum chemistry, no engine-specific input/output handling --
that stays in each engine's own SEAMM step (``orca_step`` first). See
``docs/developer_guide/campaigns/2026-08-03/bsse_architecture.rst`` for the
full design.
"""

from .fragment import Fragment  # noqa: F401
from .job_specs import (  # noqa: F401
    CLUSTER,
    FRAGMENT_ALONE,
    FRAGMENT_IN_CLUSTER,
    JobSpec,
    generate_job_specs,
)
from .combine import CPResult, JobResult, combine  # noqa: F401
from .validate import validate_fragments  # noqa: F401
from ._version import __version__  # noqa: F401

__all__ = [
    "Fragment",
    "JobSpec",
    "generate_job_specs",
    "CLUSTER",
    "FRAGMENT_IN_CLUSTER",
    "FRAGMENT_ALONE",
    "JobResult",
    "CPResult",
    "combine",
    "validate_fragments",
    "__version__",
]
