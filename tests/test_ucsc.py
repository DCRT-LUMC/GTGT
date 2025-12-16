from typing import Any

import pytest

from gtgt import Bed
from gtgt.range import Range
from gtgt.ucsc import _blocks_overlap, _track_to_range, _tracks_to_bed

payload = dict[str, Any]


@pytest.mark.parametrize(
    "tracks, expected",
    [
        # A single track
        (
            [
                {
                    "chrom": "chr11",
                    "chromStart": 10,
                    "chromEnd": 20,
                    "name": "name",
                    "blockCount": 1,
                },
            ],
            [Bed("chr11", 10, 20, name="name")],
        ),
        (
            # Two tracks with the same name
            [
                {
                    "chrom": "chr11",
                    "chromStart": 10,
                    "chromEnd": 20,
                    "name": "name",
                    "blockCount": 1,
                },
                {
                    "chrom": "chr11",
                    "chromStart": 30,
                    "chromEnd": 40,
                    "name": "name",
                    "blockCount": 1,
                },
            ],
            [
                Bed(
                    "chr11",
                    10,
                    40,
                    name="name",
                    blockCount=2,
                    blockSizes=[10, 10],
                    blockStarts=[0, 20],
                )
            ],
        ),
        (
            # Two tracks with different names do not get merged
            [
                {
                    "chrom": "chr11",
                    "chromStart": 10,
                    "chromEnd": 20,
                    "name": "name",
                    "blockCount": 1,
                },
                {
                    "chrom": "chr11",
                    "chromStart": 30,
                    "chromEnd": 40,
                    "name": "NAME",
                    "blockCount": 1,
                },
            ],
            [
                Bed("chr11", 10, 20, name="name"),
                Bed("chr11", 30, 40, name="NAME"),
            ],
        ),
        (
            # Two tracks with the same name, one with multiple blocks
            [
                # [(5, 10), (15, 25)]
                {
                    "chrom": "chr11",
                    "chromStart": 5,
                    "chromEnd": 20,
                    "name": "name",
                    "blockCount": 2,
                    "chromStarts": "0,10",  # Relative to chromStart
                    "blockSizes": "5, 10",
                },
                # [(30, 40)]
                {
                    "chrom": "chr11",
                    "chromStart": 30,
                    "chromEnd": 40,
                    "name": "name",
                    "blockCount": 1,
                },
            ],
            [
                Bed(
                    "chr11",
                    5,
                    40,
                    name="name",
                    blockCount=3,
                    blockSizes=[5, 10, 10],
                    blockStarts=[0, 10, 25],
                )
            ],
        ),
    ],
)
def test_tracks_to_bed(tracks: list[payload], expected: list[Bed]) -> None:
    # Update constant values for the tracks
    constant = {
        "reserved": "100,100,0",
        "score": 1000,
        "strand": "+",
    }
    for track in tracks:
        track.update(constant)

    # Update the constant values for the BED records
    for bed in expected:
        bed.score = 1000
        bed.strand = "+"
        bed.itemRgb = (100, 100, 0)

    assert _tracks_to_bed(tracks) == expected


@pytest.mark.parametrize(
    "track, ranges",
    [
        (
            # A single block
            {
                "chromStart": 10,
                "chromEnd": 20,
                "blockCount": 1,
                "chromStarts": "0",
                "blockSizes": "10",
            },
            [(10, 20)],
        ),
        (
            # A track with two blocks
            {
                "chromStart": 10,
                "chromEnd": 40,
                "blockCount": 2,
                "chromStarts": "0,20",  # Relative to chromStart
                "blockSizes": "10, 10",
            },
            [(10, 20), (30, 40)],
        ),
        (
            # A track with two blocks, with trailing commas
            {
                "chromStart": 10,
                "chromEnd": 40,
                "blockCount": 2,
                "chromStarts": "0,20,",  # Relative to chromStart
                "blockSizes": "10, 10,",
            },
            [(10, 20), (30, 40)],
        ),
    ],
)
def test_track_to_range(track: payload, ranges: list[Range]) -> None:
    """Test converting a single track to a list of ranges"""
    assert _track_to_range(track) == ranges


@pytest.mark.parametrize(
    "blocks, expected",
    [
        ([], False),
        ([(0, 10)], False),
        ([(0, 10), (0, 10)], True),
        ([(0, 10), (10, 20)], False),
        ([(0, 10), (10, 20), (8, 12), (0, 10)], True),
    ],
)
def test_blocks_overlap(blocks: list[Range], expected: bool) -> None:
    assert _blocks_overlap(blocks) == expected
