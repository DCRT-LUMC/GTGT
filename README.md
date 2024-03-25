[![Continous integration](https://github.com/DCRT-LUMC/GTGT/actions/workflows/ci.yml/badge.svg)](https://github.com/DCRT-LUMC/GTGT/actions/workflows/ci.yml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Checked with mypy](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![Documentation Status](https://readthedocs.org/projects/gtgt/badge/?version=latest)](https://gtgt.readthedocs.io/en/latest/?badge=latest)

# Python project
------------------------------------------------------------------------

## Documentation
The documentation is available on [http://gtgt.readthedocs.io/](http://gtgt.readthedocs.io/).

## Caching
To speed up the tool, you can use caching by either specifying a folder using `--cachedir`, or by setting the `GTGT_CACHE` environment variable.

## Human
gtgt --cachedir cache transcript ENST00000241453.12 | jq .

## Rat
gtgt --cachedir cache transcript ENSRNOT00000088957.2 | jq .
