"""The test suite, a subpackage so its module names never collide with another checkout's tests.

Run from the repository root: python -m unittest discover -s indexed_gam_pipeline_v3/tests -t .
Another pipeline package (e.g. a frozen run source) must not be loaded in the same process: only one
_fastdecode module can load per process. The goldens run every build in a subprocess (tests/golden.py).
"""
import sys

_PACKAGE = __name__.split(".")[0]
_OTHERS = sorted(name for name in sys.modules
                 if name.startswith("indexed_gam_pipeline") and name.split(".")[0] != _PACKAGE)
if _OTHERS:
    raise ImportError(f"{__name__} refuses to load next to another pipeline package ({_OTHERS[0]} is imported); "
                      "run each package's suite in its own interpreter")
