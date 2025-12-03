[![PyPi](https://img.shields.io/pypi/v/gtgt.svg)](https://pypi.org/project/gtgt)
[![Release date](https://img.shields.io/github/release-date/DCRT-LUMC/GTGT.svg)](https://github.com/DCRT-LUMC/GTGT/releases)
[![Last commit](https://img.shields.io/github/last-commit/DCRT-LUMC/GTGT.svg)](https://github.com/DCRT-LUMC/GTGT/graphs/commit-activity)
[![Tests](https://github.com/DCRT-LUMC/GTGT/actions/workflows/ci.yml/badge.svg)](https://github.com/DCRT-LUMC/GTGT/actions/workflows/ci.yml)
[![Docs](https://readthedocs.org/projects/gtgt/badge/?version=latest)](https://gtgt.readthedocs.io/en/latest/?badge=latest)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Checked with mypy](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![License](https://img.shields.io/github/license/DCRT-LUMC/GTGT.svg)](https://raw.githubusercontent.com/DCRT-LUMC/GTGT/main/LICENSE.md)

# Genetic Therapy Generator Toolkit
------------------------------------------------------------------------
## GTGT website
Try out GTGT via the website of the Dutch Center for RNA Therapeutics:
[https://gtgt.rnatherapy.nl](https://gtgt.rnatherapy.nl).

## Documentation
The documentation is available on [http://gtgt.readthedocs.io/](http://gtgt.readthedocs.io/).

## Caching
To speed up the tool, you can use caching by either specifying a folder using `--cachedir`, or by setting the `GTGT_CACHE` environment variable. This currently does not work with `gtgt api_server`.

## Human
gtgt --cachedir cache transcript ENST00000241453.12 | jq .


## Variant Information
gtgt links "NM_000094.4:c.5299G>C"

## Analyze exon skips
Use the `analyze` entry point to generate all exon skips, and compare them to the wildtype and patient transcripts

```bash
# A frameshift deletion in exon 2
$ gtgt analyze ENST00000375549.8:c.100del
{
 "wildtype": 1.0,
 "ENST00000375549.8:c.50_172del": 0.8048779330148124,
 "patient": 0.3513663563397666,
 "ENST00000375549.8:c.167_317del": 0.3513663563397666
}

# An in-frame deletion in exon 2, notice how non of the exon skips have a
# higher score than the patient
$ gtgt analyze ENST00000375549.8:c.100_102del
{
 "wildtype": 1.0,
 "patient": 0.9970458173607387,
 "ENST00000375549.8:c.50_172del": 0.8048779330148124,
 "ENST00000375549.8:c.167_317del": 0.3484121737005053
}


# An in-frame deletion that creates a STOP codon is recognized as
# as highly detrimental
$ gtgt analyze ENST00000452863.10:c.87_89del
{
 "wildtype": 1.0,
 "patient": 0.18847136926335686,
 "ENST00000452863.10:c.659_787del": 0.18847136926335686,
 "ENST00000452863.10:c.782_890del": 0.18847136926335686,
 "ENST00000452863.10:c.885_968del": 0.18847136926335686,
 "ENST00000452863.10:c.963_1019del": 0.18847136926335686,
 "ENST00000452863.10:c.1014_1116del": 0.18847136926335686,
 "ENST00000452863.10:c.1111_1267del": 0.18847136926335686,
 "ENST00000452863.10:c.1262_1357del": 0.18847136926335686,
 "ENST00000452863.10:c.1352_1450del": 0.18847136926335686
}

```

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
