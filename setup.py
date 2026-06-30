from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "fast_writer",
        ["cpp/fast_writer.cpp"],
        extra_compile_args=["-O3", "-DNDEBUG", "-std=c++17"],
    ),
]

setup(
    name="pansoma-fast-writer",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)

