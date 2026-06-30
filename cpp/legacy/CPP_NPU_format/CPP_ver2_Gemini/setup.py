# file: setup.py
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension("fast_writer", ["fast_writer.cpp"]),
]

setup(
    name="fast_writer",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)