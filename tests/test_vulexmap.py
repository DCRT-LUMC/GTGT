from typing import Any
import pytest

from gtgt.vulexmap import VulExMap, VulExMapExon, _vulexmap_key, vulexmap_description
from gtgt.mutalyzer import init_description
from mutalyzer.description import Description


def SDHD_description() -> Description:
    """SDHD, on the forward strand"""
    return init_description("ENST00000375549.8:c.=")


def WT1_description() -> Description:
    """WT1, on the reverse strand"""
    return init_description("ENST00000452863.10:c.=")


@pytest.mark.parametrize(
    "d, exon, expected",
    [
        (
            SDHD_description(),
            # Exon 2 of SDHD
            (984, 1101),
            # Key
            # chrom          start      end        strand
            ("NC_000011.10", 112087856, 112087973, "+"),
        ),
        (
            WT1_description(),
            # Exon 2 of WT1 (on the reverse strand)
            (40722, 40845),
            # Key
            # chrom          start      end        strand
            ("NC_000011.10", 32428496, 32428619, "-"),
        ),
    ],
)
def test_vulexmap_key(
    d: Description, exon: tuple[int, int], expected: tuple[Any, ...]
) -> None:
    assert _vulexmap_key(d, exon) == expected


@pytest.mark.parametrize(
    "exon_class, low_read_count, name, expected",
    [
        (
            "vulnerable",
            False,
            "Exon 1",
            "Exon 1 is a vulnerable exon (median PSI=77.8, sd=6.6).",
        ),
        (
            "resilient",
            False,
            "Exon 2",
            "Exon 2 is a resilient exon (median PSI=77.8, sd=6.6).",
        ),
        (
            "resilient",
            True,
            "Exon 3",
            "Exon 3 appears to be a resilient exon (median PSI=77.8, sd=6.6), but this is based on limited data.",
        ),
    ],
)
def test_vulexmap_description(
    exon_class: str, low_read_count: bool, name: str, expected: str
) -> None:
    v = VulExMapExon(
        median_PSI=77.77777777,
        sd_PSI=6.63333,
        exon_class=exon_class,
        low_read_count=low_read_count,
    )
    assert vulexmap_description(v, name).startswith(expected)
