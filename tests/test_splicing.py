#!/usr/bin/env python3

import pytest

from gtgt.splicing import ExonTranscript, JunctionTranscript

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
# fmt: on


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
    "novel_splice, expected_junctions",
    [
        (
            # Junction before the transcript start
            (80, 90),
            [(125, 160), (172, 208), (241, 300)],
        ),
        (
            # Junction after the transcript end
            (350, 400),
            [(125, 160), (172, 208), (241, 300)],
        ),
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
        # (
        #     # Junction from exon 3 to exon 3
        #     (165, 170),
        #     [(125, 160), (165, 170), (172, 208), (241, 300)]
        # ),
    ],
)
def test_alternative_splice_event(
    novel_splice: tuple[int, int], expected_junctions: list[tuple[int, int]]
) -> None:
    t = JunctionTranscript(
        start=100, end=348, junctions=[(125, 160), (172, 208), (241, 300)]
    )
    assert t._alternative_junctions(novel_splice) == expected_junctions
