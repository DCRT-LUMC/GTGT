#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

from glob import glob
from os.path import basename
from os.path import splitext
from pathlib import Path

from setuptools import find_packages
from setuptools import setup


this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setup(
    name='gtgt',
    version='0.2.8',
    license='AGPL-3.0',
    description='Genetic Therapy Generator Toolkit',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Redmar van den Berg',
    author_email='RedmarvandenBerg@lumc.nl',
    url='https://github.com/DCRT-LUMC/GTGT',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14',
        # uncomment if you test on these interpreters:
        # 'Programming Language :: Python :: Implementation :: IronPython',
        # 'Programming Language :: Python :: Implementation :: Jython',
        # 'Programming Language :: Python :: Implementation :: Stackless',
        'Topic :: Utilities',
    ],
    project_urls={
        'Changelog': 'https://github.com/DCRT-LUMC/GTGT/blob/main/CHANGELOG.rst',
        'Issue Tracker': 'https://github.com/DCRT-LUMC/GTGT/issues',
    },
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    python_requires='>=3.10',
    install_requires=[
        "pydantic",
        "setuptools",
        "mutalyzer>=3.1.1",
        "lark==1.2.2",
        "exonviz>=0.2.16",
    ],
    extras_require={
        "api_server": ["uvicorn", "fastapi", "mutalyzer-hgvs-parser", "httpx"],
        "webserver": ["flask"],
    },
    setup_requires=[
        'pytest-runner',
    ],
    entry_points={
        'console_scripts': [
            'gtgt=gtgt.cli:main',
        ]
    },
)
