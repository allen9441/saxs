from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "saxs_core",
        ["saxs_core.cpp"],
        cxx_std=11,
    ),
]

setup(
    name="saxs_core",
    version="0.0.1",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
