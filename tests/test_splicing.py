#!/usr/bin/env python3

from gtgt.splicing import Exon, Junction, exons_to_junctions, junctions_to_exons
import pytest

"""
These are tests for combining alternative splice events with the splice
junctions as defined in a given transcript

All tests use the following transcript definition of a transcript on the
reverse strand which starts on position 110 of a hypothetical genome

100    110      125       160    172        208    241        300    348 
v       v        v         v      v          v      v          v      v    
--------||||||||||---------||||||||----------||||||||----------||||||||
          exon4    intron3  exon3   intron2   exon2   intron1    exon1
        ^        ^         ^      ^          ^      ^          ^      ^
        r.112    r.97      r.96   r.84      r.83    r.50       r.49   r.1

"""

# fmt: off
exons = [
    (110, 125),
    (160, 172),
    (208, 241),
    (300, 348)
]
junctions = [
    (125, 160),
    (172, 208),
    (241, 300)
]

# Location of a the novel splice site in each exon/intron
splice_site_locations = {
    "exon 4": 115,
    "intron 3": 140,
    "exon 3": 165,
    "intron 2": 180,
    "exon 2": 225,
    "intron 1": 250,
    "exon 1": 333
}

# All possible alternative splice junctions that can occur in the example
# transcript
alternative_splice_junctions = [
    # Junction from exon 4 to exon 4
    (115, 120),
    # Junction from exon 4 to intron 3
    (115, 140),
    # Junction from exon 4 to exon 3
    (115, 165),
    # Junction from exon 4 to intron 2
    (115, 180),
    # Junction from exon 4 to exon 2
    (115, 225),
    # Junction from exon 4 to intron 1
    (115, 250),
    # Junction from exon 4 to exon 1
    (115, 333),
    # Junction from intron 3 to intron 3
    (140, 150),
    # Junction from intron 3 to exon 3
    (140, 165),
    # Junction from intron 3 to intron 2
    (140, 180),
    # Junction from intron 3 to exon 2
    (140, 225),
    # Junction from intron 3 to intron 1
    (140, 250),
    # Junction from intron 3 to exon 1
    (140, 333),
    # Junction from exon 3 to exon 3
    (165, 170),
    # Junction from exon 3 to intron 2
    (165, 180),
    # Junction from exon 3 to exon 2
    (165, 225),
    # Junction from exon 3 to intron 1
    (165, 250),
    # Junction from exon 3 to exon 1
    (165, 333),
    # Junction from intron 2 to intron 2
    (180, 200),
    # Junction from intron 2 to exon 2
    (180, 225),
    # Junction from intron 2 to intron 1
    (180, 250),
    # Junction from intron 2 to exon 1
    (180, 333),
    # Junction from exon 2 to exon 2
    (225, 235),
    # Junction from exon 2 to intron 1
    (225, 250),
    # Junction from exon 2 to exon 1
    (225, 333),
    # Junction from intron 1 to intron 1
    (250, 275),
    # Junction from intron 1 to exon 1
    (250, 333)
]


# Tests cases for converting between exons and junctions
exon_junction_conversion = [
        ( # Single exon
            [(10, 20)],
            10,
            20,
            []),
        ( # Two exons
            [(10, 20), (40, 50)],
            10,
            50,
            [(20, 40)]
        ),
        ( # Three exons
            [(11, 20), (40, 50), (100, 113)],
            11,
            113,
            [(20, 40), (50, 100)]
        ),
    ]
@pytest.mark.parametrize( "exons, start, end, junctions", exon_junction_conversion)
def test_exons_to_junctions(
    exons: list[Exon], start: int, end: int, junctions: list[Junction]
) -> None:
    """Test converting exons to junctions"""
    assert exons_to_junctions(exons) == (start, end, junctions)

@pytest.mark.parametrize(
    "exons, start, end, junctions", exon_junction_conversion
)
def test_junctions_to_exons(
    exons: list[Exon], start: int, end: int, junctions: list[Junction]
) -> None:
    """Test converting junctions to exons"""
    assert junctions_to_exons(start, end, junctions) == exons

# fmt: on
