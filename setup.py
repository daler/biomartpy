import os
from setuptools import setup, find_packages
import sys


long_description = """
Simple interface for accessing biomaRt from Python (via rpy2 and pandas)
"""

setup(
        name="biomartpy",
        version='0.1.1',
        install_requires=['pandas', 'rpy2'],
        packages=find_packages(),
        author="Ryan Dale",
        description=long_description,
        long_description=long_description,
        url="none",
        package_dir = {"biomartpy": "biomartpy"},
        author_email="dalerr@niddk.nih.gov",
        classifiers=['Development Status :: 4 - Beta'],
    )
