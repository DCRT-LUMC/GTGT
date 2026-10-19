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
[https://gtgt.rnatherapy.nl](https://gtgt.rnatherapy.nl). For information on
how to use GTGT, please see the online documentation on
[http://gtgt.readthedocs.io/](http://gtgt.readthedocs.io/).

## Quick start guide
The following is a quick guide on how to use GTGT locally. To speed up the
analysis of multiple variants, please see the [documentation]() on how to set
up caching for GTGT.

### Command line version
```bash

# Install GTGT
pip install gtgt

# Print the transcript information (in bed format)
gtgt transcript ENST00000375549.8

# Example frameshift variant in exon 2 of SDHD to analyze
dna=ENST00000375549.8:c.100del
rna=ENST00000375549.8:r.100del

# Analyze the variants, which will return a list of possible therapies and how
# they compare with the wildtype transcript
gtgt analyze --dna ${dna} --rna ${rna}
[
 ...
 {
  "therapy": {
   "name": "Skip exon 2",
   "hgvsc": "ENST00000375549.8:c.100del",
   "description": "The annotations based on the supplied variants, in combination with skipping exon 2.\n",
   "variants": [],
   "figure": null,
   "hgvsr": "ENST00000375549.8:r.53_169del",
   "hgvsp": "ENST00000375549.8(ENSP00000364699):p.(Leu19_Ser57del)"
  },
  "comparison": [
   {
    "name": "Exons",
    "percentage": 0.912621359223301,
    "basepairs": "1222/1339"
   },
   {
    "name": "Coding exons",
    "percentage": 0.7520833333333333,
    "basepairs": "361/480"
   }
  ]
 },
 ...
]

# To see the exact features of the transcript after skipping exon 2 (see above), use
gtgt export --dna ${dna} --rna "ENST00000375549.8:r.53_169del" --protein
```

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
