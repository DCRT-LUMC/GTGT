#!/usr/bin/env python3

import pytest

from gtgt.splicing import (
    ExonTranscript,
    Junction,
    JunctionTranscript,
    _must_be_increasing,
)

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



# Tests cases for converting between exons and junctions
exon_junction_conversion = [
    (  # Single exon
        ExonTranscript([(10, 20)]),
        JunctionTranscript(10, 20, []),
    ),
    (  # Two exons
        ExonTranscript([(10, 20), (40, 50)]),
        JunctionTranscript(10, 50, [(20, 40)]),
    ),
    (  # Three exons
        ExonTranscript([(11, 20), (40, 50), (100, 113)]),
        JunctionTranscript(11, 113, [(20, 40), (50, 100)]),
    ),
]


@pytest.mark.parametrize("exons, junctions", exon_junction_conversion)
def test_exons_to_junctions(
    exons: ExonTranscript, junctions: JunctionTranscript
) -> None:
    """Test converting exons to junctions"""
    assert exons.to_junctions() == junctions


@pytest.mark.parametrize("exons, junctions", exon_junction_conversion)
def test_junctions_to_exons(
    exons: ExonTranscript, junctions: JunctionTranscript
) -> None:
    """Test converting junctions to exons"""
    assert junctions.to_exons() == exons


@pytest.mark.parametrize(
    "novel_junction, expected_junctions",
    [
        (
            # Junction before the first regular junction
            (115, 120),
            [(115, 120), (125, 160), (172, 208), (241, 300)],
        ),
        (
            # Junction after the last regular junction
            (320, 340),
            [(125, 160), (172, 208), (241, 300), (320, 340)],
        ),
        (
            # Junction fully in exon 4
            (115, 120),
            [(115, 120), (125, 160), (172, 208), (241, 300)],
        ),
        (
            # Junction from exon 4 to intron 3
            (115, 140),
            [(115, 140), (172, 208), (241, 300)],
        ),
        (
            # Junction from exon 4 to exon 3
            (115, 165),
            [(115, 165), (172, 208), (241, 300)],
        ),
        (
            # Junction from exon 4 to intron 2
            (115, 180),
            [(115, 180), (241, 300)],
        ),
        (
            # Junction from exon 4 to exon 2
            (115, 225),
            [(115, 225), (241, 300)],
        ),
        (
            # Junction from exon 4 to intron 1
            (115, 250),
            [(115, 250)],
        ),
        (
            # Junction from exon 4 to exon 1
            (115, 333),
            [(115, 333)],
        ),
        (
            # Junction from intron 3 to intron 3
            (140, 150),
            [(140, 150), (172, 208), (241, 300)],
        ),
        (
            # Junction from intron 3 to exon 3
            (140, 165),
            [(140, 165), (172, 208), (241, 300)],
        ),
        (
            # Junction from intron 3 to intron 2
            (140, 180),
            [(140, 180), (241, 300)],
        ),
        (
            # Junction from intron 3 to exon 2
            (140, 225),
            [(140, 225), (241, 300)],
        ),
        (
            # Junction from intron 3 to intron 1
            (140, 250),
            [(140, 250)],
        ),
        (
            # Junction from intron 3 to exon 1
            (140, 333),
            [(140, 333)],
        ),
        (
            # Junction from exon 3 to exon 3
            (165, 170),
            [(125, 160), (165, 170), (172, 208), (241, 300)],
        ),
        (
            # Junction from exon 3 to intron 2
            (165, 180),
            [(125, 160), (165, 180), (241, 300)],
        ),
        (
            # Junction from exon 3 to exon 2
            (165, 225),
            [(125, 160), (165, 225), (241, 300)],
        ),
        (
            # Junction from exon 3 to intron 1
            (165, 250),
            [(125, 160), (165, 250)],
        ),
        (
            # Junction from exon 3 to exon 1
            (165, 333),
            [(125, 160), (165, 333)],
        ),
        (
            # Junction from intron 2 to intron 2
            (180, 200),
            [(125, 160), (180, 200), (241, 300)],
        ),
        (
            # Junction from intron 2 to exon 2
            (180, 225),
            [(125, 160), (180, 225), (241, 300)],
        ),
        (
            # Junction from intron 2 to intron 1
            (180, 250),
            [(125, 160), (180, 250)],
        ),
        (
            # Junction from intron 2 to exon 1
            (180, 333),
            [(125, 160), (180, 333)],
        ),
        (
            # Junction from exon 2 to exon 2
            (225, 235),
            [(125, 160), (172, 208), (225, 235), (241, 300)],
        ),
        (
            # Junction from exon 2 to intron 1
            (225, 250),
            [(125, 160), (172, 208), (225, 250)],
        ),
        (
            # Junction from exon 2 to exon 1
            (225, 333),
            [(125, 160), (172, 208), (225, 333)],
        ),
        (
            # Junction from intron 1 to intron 1
            (250, 275),
            [(125, 160), (172, 208), (250, 275)],
        ),
        (
            # Junction from intron 1 to exon 1
            (250, 333),
            [(125, 160), (172, 208), (250, 333)],
        ),
        (
            # Skip exon 3
            (125, 208),
            [(125, 208), (241, 300)],
        ),
        (
            # Skip exon 2
            (172, 300),
            [(125, 160), (172, 300)],
        ),
    ],
)
def test_alternative_splice_event(
    novel_junction: Junction, expected_junctions: list[Junction]
) -> None:
    t = JunctionTranscript(
        start=100, end=348, junctions=[(125, 160), (172, 208), (241, 300)]
    )
    assert t._alternative_junctions(novel_junction) == expected_junctions

@pytest.mark.parametrize(
    "novel_junction",
    [
        # Junction before the transcript start
        (80, 90),
        # Junction after the transcript end
        (350, 400),
        # Junction from before the start of the transcript
        (325, 400),
        # Junction to beyond the end of the transcript
        (90, 115)
    ],
)
def test_alternative_splice_event_ignored_junctions( novel_junction: Junction) -> None:
    """ These alternative junctions are outside the exons and should be ignored"""
    junctions = [(125, 160), (172, 208), (241, 300)]
    t = JunctionTranscript(
        start=100, end=348, junctions=junctions
    )
    assert t._alternative_junctions(novel_junction) == junctions

def test_junctions_in_order() -> None:
    """ Raise an error when the junction end is not after the start"""
    with pytest.raises(ValueError):
        t = JunctionTranscript(100, 200, junctions=[])
        t._alternative_junctions((10, 10))


@pytest.mark.parametrize(
    "junctions",
    [
        # The end is not after the start
        [(0,0)],
        # The junctions are not in order
        [(10, 20), (4, 8)],
        # The junctions overlap
        [(10, 20), (19, 21)],
        # The junctions overlap and are out of order
        [(19, 21), (10, 20)],
    ]
)
def test_invalid_junctions(junctions: list[Junction]) -> None:
    """ Check that we raise an error when we find invalid junctions

    A list of junctions is invalid if:
        - The the end is not after the start for any junction
        - The junctions are not in order
        - The junctions overlap
    """
    with pytest.raises(ValueError):
        _must_be_increasing(junctions)
