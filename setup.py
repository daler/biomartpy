import ez_setup
ez_setup.use_setuptools()

import os
import sys
from setuptools import setup


long_description = """
Simple interface for accessing biomaRt from Python (via rpy2 and pandas)
"""

setup(
        name="biomartpy",
        version=0.1,
        install_requires=['pandas', 'rpy2'],
        packages=['biomartpy',
                  ],
        author="Ryan Dale",
        description=long_description,
        long_description=long_description,
        url="none",
        package_dir = {"biomartpy": "biomartpy"},
        author_email="dalerr@niddk.nih.gov",
        classifiers=['Development Status :: 4 - Beta'],
    )
