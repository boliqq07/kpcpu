from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages

ext_modules = []

# pybind11
ext_modules.append(Pybind11Extension("kpcpu", ["hskpt/kpcpu.cpp"], language="c++"))
include_dirs = [r'/home/iap13/eigen-3.3.9']


def main():
    setup(name="kpcpu",
          version="0.0.5",
          cmdclass={"build_ext": build_ext},
          description="Python / C library ",
          author="boliqq07",
          author_email="98988989@qq.com",
          packages=find_packages(exclude=[".tests", ".tests.", "tests.", "tests"]),
          platforms=[
              "Windows",
              "Unix",
          ],
          classifiers=[
              "Development Status :: 4 - Beta",
              "Intended Audience :: Science/Research",
              "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
              "Natural Language :: English",
              "Operating System :: Microsoft :: Windows",
              "Operating System :: Unix",
              "Programming Language :: Python :: 3.6",
              "Programming Language :: Python :: 3.7",
              "Programming Language :: Python :: 3.8",
              "Programming Language :: Python :: 3.9",
          ],
          include_dirs=include_dirs,
          ext_modules=ext_modules)


# python setup.py bdist_wheel


if __name__ == "__main__":
    main()
