#!/usr/bin/env python3

import logging
from dataclasses import dataclass
from typing import Any

from mutalyzer.description import Description

from gtgt.mutalyzer import chrom_to_nc, get_chrom_name, get_offset

logger = logging.getLogger(__name__)


@dataclass
class VulExMapExon:
    median_PSI: float
    sd_PSI: float
    exon_class: str
    low_read_count: bool


VulExDict = dict[Any, tuple[float, float, str, bool]]


class VulExMap:
    """A caching object to access VulExMap data"""

    def __init__(self, path: str | None = None):
        self.path = path
        self._data: VulExDict | None = None

    def __getitem__(self, key: Any) -> VulExMapExon | None:
        # If path is not defined, we always return None
        if not self.path:
            return None

        # Read the file only once
        if self.path is not None and self._data is None:
            try:
                self._data = self.read_vulex_file(self.path)
            except:
                logger.warning(f"Failed to read data for {self}")
                self._data = dict()

        values = self._data.get(key)
        if values:
            return VulExMapExon(*values)
        else:
            return None

    def __repr__(self) -> str:
        return f"VulExMap(path={self.path})"

    def read_vulex_file(self, fname: str) -> VulExDict:
        results = dict()
        with open(fname) as fin:
            header = next(fin).strip().split("\t")

            for line in fin:
                spline = line.strip().split("\t")
                d = {k: v for k, v in zip(header, spline)}

                # Prepare the lookup key
                chr = chrom_to_nc(d.pop("chr"))
                start = int(d.pop("start")) - 1
                end = int(d.pop("end"))
                strand = d.pop("strand")

                key = (chr, start, end, strand)

                # Convert values to the correct type
                median_PSI = float(d["medianPSI"])
                sd_PSI = float(d["sd_PSI"])
                exon_class = d["class"]
                low_read_count = d["Mean_read_count_ano"] == "*"

                results[key] = (median_PSI, sd_PSI, exon_class, low_read_count)
        return results


def _vulexmap_key(d: Description, exon: tuple[int, int]) -> tuple[Any, ...]:
    """Determine the VulExMap key from a description and exon coordinates"""
    # Create the key to look up VulExMap data
    chr = get_chrom_name(d)
    offset = get_offset(d)
    start = exon[0] + offset
    end = exon[1] + offset
    strand = "-" if d.is_inverted() else "+"
    return (chr, start, end, strand)


def lookup_vulexmap(
    d: Description, exon: tuple[int, int], V: VulExMap = VulExMap()
) -> VulExMapExon | None:
    key = _vulexmap_key(d, exon)

    return V[key]


def vulexmap_description(v: VulExMapExon, name: str) -> str:
    """Create a human readable description based on the VulExMap data"""
    # fmt: off
    if v.low_read_count:
        return (
            f"{name} appears to be a {v.exon_class} exon "
            f"(median PSI={v.median_PSI:.1f}, "
            f"sd={v.sd_PSI:.1f}), but this is based on limited data."
        )
    return (
            f"{name} is a {v.exon_class} exon "
            f"(median PSI={v.median_PSI:.1f}, "
            f"sd={v.sd_PSI:.1f})."
    )
    # fmt: on


def main(vulex: str) -> None:

    V = VulExMap(path=vulex)
    # V = VulEx()
    for (
        hgvs
    ) in "ENST00000375549.8:c.100del NC_000023.11(NM_004006.3):c.2500del".split():
        d = init_description(hgvs)
        print(d)
        offset = get_offset(d)
        chr = get_chrom_name(d)
        strand = "-" if d.is_inverted() else "+"

        for exon in get_exons(d, in_transcript_order=True):
            start, end = exon
            start += offset
            end += offset

            key = (chr, start, end, strand)

            print(hgvs, key, V[key])


if __name__ == "__main__":
    import sys

    from gtgt.mutalyzer import get_chrom_name, get_exons, get_offset, init_description

    main(sys.argv[1])
