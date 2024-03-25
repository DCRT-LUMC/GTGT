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


# Disclaimer
Copyright© 2023 LUMC (https://www.lumc.nl)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

By accessing and using the program in any manner (including copying, modifying
or redistributing the program), you accept and agree to the applicability of
the GNU Affero General Public License. You can find and read this license on
GNU Affero General Public License - GNU Project - Free Software Foundation.

In case of questions, you can contact us at DCRT@LUMC.nl.
