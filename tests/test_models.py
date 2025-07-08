from gtgt.models import BedModel, TranscriptId, TranscriptModel, TranscriptModel
from gtgt.mutalyzer import HGVS
from gtgt.bed import Bed
from gtgt.transcript import Transcript
import pytest
from typing import Dict, Union, Tuple
from pydantic import ValidationError

payload = Dict[str, Union[str, int]]


@pytest.fixture
def ucsc() -> payload:
    return {
        "chrom": "chr1",
        "chromStart": 1000,
        "chromEnd": 2000,
        "name": "ENST00000.12",
        "score": 0,
        "strand": "-",
        "thickStart": 1000,
        "thickEnd": 2000,
        "blockCount": 2,
        "blockSizes": "200,700,",
        "chromStarts": "0,300,",
        "random_field": "some nonsense",
    }


def test_model_from_ucsc(ucsc: payload) -> None:
    """Test creating a BedModel from UCSC payload"""
    bm = BedModel.from_ucsc(ucsc)

    expected = Bed(
        "chr1",
        1000,
        2000,
        "ENST00000.12",
        0,
        "-",
        thickStart=1000,
        thickEnd=2000,
        blockCount=2,
        blockSizes=[200, 700],
        blockStarts=[0, 300],
    )
    new_bed = bm.to_bed()

    assert new_bed == expected


def test_Bed_from_model(ucsc: payload) -> None:
    """Test creating a Bed object from BedModel"""
    bm = BedModel.from_ucsc(ucsc)

    bed = bm.to_bed()
    assert bed.chrom == "chr1"
    assert bed.itemRgb == (0, 0, 0)
    assert bed.blockStarts == [0, 300]


def test_Bed_validation(ucsc: payload) -> None:
    # Blocks must be in ascending order
    ucsc["chromStarts"] = "300,0"
    ucsc["blockSizes"] = "700,200"
    with pytest.raises(ValueError):
        BedModel.from_ucsc(ucsc)


def test_BedModel_from_bed() -> None:
    bed = Bed("chr1", 0, 10)
    bm = BedModel.from_bed(bed)
    assert bm.chrom == "chr1"
    assert bm.blocks == [(0, 10)]


def test_transcript_model() -> None:
    bed = Bed("chr1", 0, 10)
    bm = BedModel.from_bed(bed)
    tm = TranscriptModel(exons=bm, cds=bm)
    # Test converting a TranscriptModel to a Transcript
    transcript = tm.to_transcript()
    # Test that the "coding" region has been set in the new Transcript
    assert transcript.coding == Bed("chr1", 0, 10, name="Coding exons")


def test_Transcript_from_model() -> None:
    """
    GIVEN a Transcript
    WHEN we create a TranscriptModel out of it
    THEN it must match the expected TranscriptModel
    """
    # Create the TranscriptModel
    bed = Bed("chr1", 0, 10)
    bm = BedModel.from_bed(bed)
    expected = TranscriptModel(exons=bm, cds=bm)

    # Create the Transcript we want to convert to a TranscriptModel
    bed = Bed("chr1", 0, 10)
    ts = Transcript(exons=bed, cds=bed)

    assert TranscriptModel.from_transcript(ts) == expected


def test_HGVS_model_valid() -> None:
    """
    GIVEN a valid HGVS description
    WHEN we make an HGVS object out of it
    THEN there should be no error
    """
    HGVS(description="NM_000094.4:c.5299G>C")


INVALID_HGVS = [
    "NM_000094.4:c.5299G>",
    "NM_000094.4>",
    "NM_000094",
]


@pytest.mark.parametrize("description", INVALID_HGVS)
def test_HGVS_model_invalid(description: str) -> None:
    """
    GIVEN an invalid HGVS description
    WHEN we make an HGVS object out of itemRgb
    THEN we should get a ValidationError
    """
    with pytest.raises(ValidationError):
        HGVS(description=description)


UNSUPPORTED = [
    # No variant present
    ("ENST:c.="),
    # Variant outside of the CDS
    ("ENST:c.-10del"),
    ("ENST:c.*10del"),
    # Intronic variants
    ("ENST:c.10+10del"),
    # Not a deletion
    ("ENST:c.10A>T"),
    # Different transcript ID
    ("ENST001:c.10del"),
]


@pytest.mark.parametrize("deletion", UNSUPPORTED)
def test_HGVS_model_add_deletion_invalid(deletion: str) -> None:
    """
    GIVEN a deletion to add to a HGVS variant
    WHEN the deletion is of an unsupported type
    THEN we raise an error
    """
    variant = HGVS(description="ENST:c.10A>T")

    with pytest.raises(NotImplementedError):
        variant.apply_deletion(HGVS(description=deletion))


UNSUPPORTED_DELETION = [
    # Partially overlaps
    ("ENST:c.10_14del"),
]


@pytest.mark.parametrize("indel", UNSUPPORTED_DELETION)
def test_HGVS_model_add_insertion_deletion_invalid(indel: str) -> None:
    """
    GIVEN a deletion to add to a HGVS variant which is an indel
    WHEN the deletion is of an unsupported type
    THEN we raise an error
    """
    variant = HGVS(description="ENST:c.10_15delinsATCG")

    with pytest.raises(NotImplementedError):
        variant.apply_deletion(HGVS(description=indel))


POSITIONS = [
    # Variant, positions
    ("ENST:c.10", (10, 10)),
    ("ENST:c.10A>T", (10, 10)),
    ("ENST:c.10_11insATC", (10, 11)),
    ("ENST:c.10_12del", (10, 12)),
]


@pytest.mark.parametrize("description, expected", POSITIONS)
def test_HGVS_model_positions(description: str, expected: Tuple[int, int]) -> None:
    """
    GIVEN a HGVS description
    WHEN we determine the position
    THEN we should get a range for the position
    """
    variant = HGVS(description=description)
    assert variant.position == expected


SUPPORTED = [
    # Deletion, expected
    # Deletion fully overlaps
    ("ENST:c.8_15del", "ENST:c.8_15del"),
    # Deletion fully overlaps, but starts at the variant
    ("ENST:c.10_15del", "ENST:c.10_15del"),
    # Deletion fully overlaps, but ends at the variant
    ("ENST:c.8_10del", "ENST:c.8_10del"),
    # Deletion is before the variant
    ("ENST:c.8del", "ENST:c.[8del;10A>T]"),
    # Deletion is after the variant
    ("ENST:c.11del", "ENST:c.[10A>T;11del]"),
]


@pytest.mark.parametrize("deletion, expected", SUPPORTED)
def test_HGVS_model_add_deletion(deletion: str, expected: str) -> None:
    """
    GIVEN a deletion to add to a HGVS variant
    WHEN the deletion is applied to the variant
    THEN the variant is updated to include the deletion
    """
    variant = HGVS(description="ENST:c.10A>T")
    variant.apply_deletion(HGVS(description=deletion))
    assert variant.description == expected


SMALLER_DELETION = [
    # Overlap start
    "10del",
    # Overlap end
    "20del",
    # Start till middle
    "10_15del",
    # Middle to end
    "15_20del",
    # Full overlap
    "10_20del",
]


@pytest.mark.parametrize("small_del", SMALLER_DELETION)
def test_HGVS_model_add_smaller_deletion(small_del: str) -> None:
    """
    GIVEN a deletion to add to a HGVS variant, which is itself a bigger deletion
    WHEN the deletion is applied to the variant
    THEN the variant should remain unchanged
    """
    variant = HGVS(description="ENST:c.10_20del")
    deletion = HGVS(description=f"ENST:c.{small_del}")
    variant.apply_deletion(deletion)
    assert variant.description == "ENST:c.10_20del"


DELETION = [
    # Deletion, expected
    ("9del", "ENST:c.[9del;10_11insATC]"),
    ("10del", "ENST:c.[10del;10_11insATC]"),
    ("11del", "ENST:c.[10_11insATC;11del]"),
    ("9_10del", "ENST:c.[9_10del;10_11insATC]"),
    ("11_12del", "ENST:c.[10_11insATC;11_12del]"),
    ("9_11del", "ENST:c.9_11del"),
    ("10_11del", "ENST:c.10_11del"),
    ("10_12del", "ENST:c.10_12del"),
    ("9_12del", "ENST:c.9_12del"),
]


@pytest.mark.parametrize("deletion, expected", DELETION)
def test_HGVS_model_add_variant_is_insertion(deletion: str, expected: str):
    """
    GIVE a deletion to add to an HGVS insertion
    WHEN the deletion is applied to the variant
    THEN the result should be both variants combined
    """
    variant = HGVS(description="ENST:c.10_11insATC")
    deletion = HGVS(description=f"ENST:c.{deletion}")
    variant.apply_deletion(deletion)
    assert variant.description == expected


DELETION = [
    # deletion, expected
    ("9del", "ENST:c.[9del;10_10delinsA]"),
    ("10del", "ENST:c.10del"),
    ("11del", "ENST:c.[10_10delinsA;11del]"),
    ("9_10del", "ENST:c.9_10del"),
    ("10_11del", "ENST:c.10_11del"),
    # (),
]


@pytest.mark.parametrize("deletion, expected", DELETION)
def test_HGVS_model_add_variant_is_one_bp_delins(deletion: str, expected: str):
    """
    GIVE a deletion to add to an HGVS indel
    WHEN the deletion is applied to the variant
    THEN the result should be both variants combined
    """
    variant = HGVS(description="ENST:c.10_10delinsA")
    deletion = HGVS(description=f"ENST:c.{deletion}")
    variant.apply_deletion(deletion)
    assert variant.description == expected


DELETION = [
    # deletion, expected
    ("9del", "ENST:c.[9del;10_12delinsGGT]"),
    ("10_12del", "ENST:c.10_12del"),
    ("9_12del", "ENST:c.9_12del"),
    ("10_13del", "ENST:c.10_13del"),
]


@pytest.mark.parametrize("deletion, expected", DELETION)
def test_HGVS_model_add_variant_is_delins(deletion: str, expected: str):
    """
    GIVE a deletion to add to an HGVS indel
    WHEN the deletion is applied to the variant
    THEN the result should be both variants combined
    """
    variant = HGVS(description="ENST:c.10_12delinsGGT")
    deletion = HGVS(description=f"ENST:c.{deletion}")
    variant.apply_deletion(deletion)
    assert variant.description == expected


VALID_TRANSCRIPT_ID = [
    "ENST00000296930.10",
]


@pytest.mark.parametrize("id", VALID_TRANSCRIPT_ID)
def test_TranscriptId_valid(id: str) -> None:
    TranscriptId(id=id)


INVALID_TRANSCRIPT_ID = [
    "ENST00000296930",
    "ENST00000296930.10:c.100A>T",
]


@pytest.mark.parametrize("id", INVALID_TRANSCRIPT_ID)
def test_TranscriptId_invalid(id: str) -> None:
    with pytest.raises(ValidationError):
        TranscriptId(id=id)
