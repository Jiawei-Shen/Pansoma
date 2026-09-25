"""The v3 test suite, a subpackage so its module names never collide with v2's top-level tests.

Run from the repository root: python -m unittest discover -s indexed_gam_pipeline_v3/tests -t .
Another pipeline package (v2) must not be loaded in the same process: only one _fastdecode
module can load per process. The goldens run v2 only in subprocesses (tests/golden.py).
"""
import sys

_PACKAGE = __name__.split(".")[0]
_OTHERS = sorted(name for name in sys.modules
                 if name.startswith("indexed_gam_pipeline") and name.split(".")[0] != _PACKAGE)
if _OTHERS:
    raise ImportError(f"{__name__} refuses to load next to another pipeline package ({_OTHERS[0]} is imported); "
                      "run the v2 and v3 suites in separate interpreters")
