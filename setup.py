#!/usr/bin/env python

import io
import glob
import os

import setuptools


def read(*names, **kwargs):
    return io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8"),
    ).read()


setuptools.setup(
    name="bcl2fastr-py",
    version="0.1.0",
    license="MIT License",
    description="Python scripts for generating CBCL files for testing",
    long_description=read("README.md"),
    author="James Webber",
    author_email="james.webber@czbiohub.org",
    url="https://github.com/czbiohub/bcl2fastr",
    packages=setuptools.find_packages("pysrc"),
    package_dir={"": "pysrc"},
    py_modules=[
        os.path.splitext(os.path.basename(path))[0] for path in glob.glob("pysrc/*.py")
    ],
    include_package_data=True,
    zip_safe=False,
    install_requires=["click >= 6.7", "numpy >= 1.15.0"],
    entry_points={
        "console_scripts": ["make_test_cbcls = bcl2fastr.make_test_cbcls:main"]
    },
)
