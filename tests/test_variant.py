import pytest
from mutalyzer.description import Description
from gtgt.mutalyzer import (
    Variant,
    CdotVariant,
    combine_variants_deletion,
    _init_model,
    mutation_to_cds_effect,
    _cdot_to_internal_delins,
    to_cdot_hgvs,
)

from pathlib import Path
from typing import List, Sequence, Tuple, Any


# Setup fixtures for mutalyzer retriever
def _retrieve_raw(
    reference_id: str,
    reference_source: Any = None,
    reference_type: Any = None,
    size_off: bool = True,
    configuration_path: Any = None,
    timeout: int = 1,
) -> Tuple[str, str, str]:
    if reference_type == "fasta":
        return _get_content("data/" + reference_id + ".fasta"), "fasta", "ncbi"
    elif reference_id.startswith("LRG_"):
        return _get_content("data/" + reference_id), "lrg", "lrg"
    else:
        return _get_content("data/" + reference_id + ".gff3"), "gff3", "ncbi"


def _get_cds_to_mrna(cds_id: Any, timeout: Any = 10) -> None:
    return None


def _get_content(relative_location: str) -> str:
    data_file = Path(__file__).parent.joinpath(relative_location)
    try:
        with open(str(data_file), "r") as file:
            content = file.read()
    except FileNotFoundError:
        raise RuntimeError({}, [])
    return content


@pytest.fixture(autouse=True)
def mock_env(monkeypatch: Any) -> None:
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", _retrieve_raw)
    monkeypatch.setattr("mutalyzer.description.get_cds_to_mrna", _get_cds_to_mrna)


@pytest.fixture()
def SDHD_description() -> Description:
    """SDHD, on the forward strand"""
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    return d


@pytest.fixture()
def WT1_description() -> Description:
    """WT1, on the reverse strand"""
    d = Description("ENST00000452863.10:c.=")
    _init_model(d)
    return d


def test_Variant_class_str() -> None:
    """Test converting a Variant to string"""
    v = Variant(10, 11, "ATG")
    assert str(v) == "Variant(start=10, end=11, inserted=ATG, deleted=)"


def test_Variant_class_str_snp() -> None:
    """SNPS are special, since they contain the inserted sequence"""
    # 10A>T
    v = Variant(10, 11, "T", "A")
    assert str(v) == "Variant(start=10, end=11, inserted=T, deleted=A)"


def test_Variant_class_to_model_positions() -> None:
    """Test converting a variant to model"""
    v = Variant(10, 11, "ATG")
    model = v.to_model()

    assert model["location"]["start"]["position"] == 10
    assert model["location"]["end"]["position"] == 11
    assert model["inserted"][0]["sequence"] == "ATG"

    # Deleted entry is missing for deletions
    assert "deleted" not in model


def test_Variant_class_to_model_snp() -> None:
    # 10 A>T
    v = Variant(10, 11, "T", "A")
    model = v.to_model()

    assert model["inserted"][0]["sequence"] == "T"
    assert model["deleted"][0]["sequence"] == "A"


def test_Variant_deleted_snp_only() -> None:
    """Test that deleted is only defined for SNPs, not larger indels

    (To match the behaviour of Mutalyzer)
    """
    with pytest.raises(ValueError):
        # 10_12delinsGG
        Variant(10, 12, "GG", "AT")


def test_Variant_class_no_inserted_sequence() -> None:
    """
    GIVEN a Variant with an empty inserted sequence
    WHEN we convert to a model
    THEN inserted should be an empty list
    """
    v = Variant(10, 11, "")
    model = v.to_model()

    assert model["inserted"] == []


def test_Variant_class_end_after_start() -> None:
    """Ensure that end is after start"""
    with pytest.raises(ValueError):
        Variant(11, 10, "")


ORDERING = [
    # Ends are touching
    (Variant(1, 3), Variant(3, 5)),
    # Gap between variants
    (Variant(0, 1), Variant(2, 4)),
]


@pytest.mark.parametrize("a, b", ORDERING)
def test_Variant_class_relative_positions(a: Variant, b: Variant) -> None:
    """Variant a is before variant b"""
    assert a.before(b)
    assert b.after(a)


INSIDE = [
    # Variants are inside themselves
    (Variant(0, 3), Variant(0, 3)),
    # Smaller variant a is inside b
    (Variant(1, 3), Variant(0, 3)),
    (Variant(0, 2), Variant(0, 3)),
    (Variant(1, 2), Variant(0, 3)),
]


@pytest.mark.parametrize("a, b", INSIDE)
def test_Variant_class_inside(a: Variant, b: Variant) -> None:
    """Variant a is inside variant b"""
    assert a.inside(b)


NOT_INSIDE = [
    # a starts outside of b
    (Variant(0, 3), Variant(1, 3)),
    # a ends outside of b
    (Variant(1, 4), Variant(1, 3)),
    # b is inside of a
    (Variant(0, 3), Variant(1, 2)),
    # a before b
    (Variant(0, 3), Variant(3, 5)),
    # a after b
    (Variant(3, 5), Variant(1, 3)),
]


@pytest.mark.parametrize("a, b", NOT_INSIDE)
def test_Variant_class_not_inside(a: Variant, b: Variant) -> None:
    """Variant a is not inside variant b"""
    assert not a.inside(b)


OVERLAP = [
    # b ends inside a
    (Variant(2, 5), Variant(1, 3)),
    # b fully inside a
    (Variant(2, 5), Variant(3, 4)),
    # b fully inside a, ends at end
    (Variant(2, 5), Variant(3, 5)),
    # b ends inside a
    (Variant(2, 5), Variant(2, 4)),
    # b starts in a, extends after
    (Variant(2, 5), Variant(3, 6)),
    # b start before a
    (Variant(2, 5), Variant(1, 5)),
    # a is inside b
    (Variant(2, 5), Variant(1, 6)),
    # a is inside b
    (Variant(2, 5), Variant(2, 6)),
    # a equals b
    (Variant(2, 5), Variant(2, 5)),
]


@pytest.mark.parametrize("a, b", OVERLAP)
def test_Variant_class_overlap(a: Variant, b: Variant) -> None:
    """Variant a and b overlap"""
    assert a.overlap(b)
    assert b.overlap(a)


NO_OVERLAP = [
    # a is after b
    (Variant(2, 5), Variant(1, 2)),
    # a is before b
    (Variant(2, 5), Variant(5, 6)),
]


@pytest.mark.parametrize("a, b", NO_OVERLAP)
def test_Variant_class_no_overlap(a: Variant, b: Variant) -> None:
    """Variant a and b do not overlap"""
    assert not a.overlap(b)
    assert not b.overlap(a)


def test_Variant_class_order() -> None:
    """Test sorting variants by start position"""
    v1 = Variant(10, 11)
    v2 = Variant(0, 1)
    assert sorted([v1, v2]) == [v2, v1]


def test_Variant_class_order_error() -> None:
    """Test error when sorting overlapping variants"""
    v1 = Variant(10, 11)
    v2 = Variant(10, 11)
    with pytest.raises(ValueError):
        sorted([v1, v2])


def test_Variant_from_model() -> None:
    """Test creating a variant from a mutalyzer delins model"""
    # fmt: off
    delins_model = {
        "location": {
            "start": {
                "position": 0
            },
            "end": {
                "position": 2
            }
        },
        "inserted": [
            {"sequence": "ATC"}
        ]
    }
    # fmt: on
    v = Variant.from_model(delins_model)
    assert v.start == 0
    assert v.end == 2
    assert v.inserted == "ATC"


def test_Variant_from_model_deletion() -> None:
    """Test creating a variant from a mutalyzer delins model

    If nothing is inserted, the inserted sequence is an emtpy list, but should
    be made an empty string in the _Variant object
    """
    # fmt: off
    delins_model = {
        "location": {
            "start": {
                "position": 0
            },
            "end": {
                "position": 2
            }
        },
        "inserted": [
            {"sequence": []}
        ]
    }
    # fmt: on
    v = Variant.from_model(delins_model)
    assert v.start == 0
    assert v.end == 2
    # Test that the emtpy sequence is a string, not a list
    assert v.inserted == ""


def test_combine_variants_deletion_empty() -> None:
    variants: List[Variant] = list()
    deletion = Variant(0, 10)
    assert combine_variants_deletion(variants, deletion) == [deletion]


# fmt: off
COMBINE = [
    ( # The variants are out of order
        # Variants
        [Variant(5, 7), Variant(2, 4)],
        # Deletion
        Variant(10, 11),
        # Expected
        [Variant(2,4), Variant(5, 7), Variant(10, 11)],
    ),
    ( # The deletion is before both variants
        # Variants
        [Variant(2, 4), Variant(5, 7)],
        # Deletion
        Variant(0, 1),
        # Expected
        [Variant(0, 1), Variant(2,4), Variant(5, 7)],
    ),
    ( # The deletion is between the variants
        # Variants
        [Variant(2, 4), Variant(5, 7)],
        # Deletion
        Variant(4, 5),
        # Expected
        [Variant(2,4), Variant(4, 5), Variant(5, 7)],
    ),
    ( # The deletion is after both variants
        # Variants
        [Variant(2, 4), Variant(5, 7)],
        # Deletion
        Variant(10, 11),
        # Expected
        [Variant(2,4), Variant(5, 7), Variant(10, 11)],
    ),
    ( # The deletion contains the first variant
        # Variants
        [Variant(2, 4), Variant(5, 7)],
        # Deletion
        Variant(1, 5),
        # Expected
        [Variant(1, 5), Variant(5, 7)],
    ),
    ( # The deletion contains the second variant
        # Variants
        [Variant(2,4), Variant(5, 7)],
        # Deletion
        Variant(4, 11),
        # Expected
        [Variant(2,4), Variant(4, 11)],
    ),
    ( # The deletion contains both variants
        # Variants
        [Variant(2,4), Variant(5, 7)],
        # Deletion
        Variant(2, 7),
        # Expected
        [Variant(2, 7)],
    ),
]
# fmt: on


@pytest.mark.parametrize("variants, deletion, expected", COMBINE)
def test_combine_variants_deletion(
    variants: Sequence[Variant], deletion: Variant, expected: Sequence[Variant]
) -> None:
    combined = combine_variants_deletion(variants, deletion)
    assert combined == expected


def test_combine_variants_deletion_variants_overlap_eachother() -> None:
    """Test that we get a value error if the variants overlap"""
    variants = [Variant(0, 2), Variant(1, 3)]
    deletion = Variant(10, 11)
    with pytest.raises(ValueError):
        combine_variants_deletion(variants, deletion)


def test_combine_variants_deletion_variants_partially_overlap_deletion() -> None:
    """Test that we get a value error if one the variants partially overlaps the deletion"""
    variants = [Variant(2, 4)]
    deletion = Variant(3, 11)
    with pytest.raises(ValueError):
        combine_variants_deletion(variants, deletion)


# fmt: off
FORWARD_MUTATIONS_VARIANT = [
    # HGVS, coordinates on the genome,
    # A simple missense that changes a single amino acids
    (
        "13T>A",
        [Variant(start=47, end=48, inserted="A", deleted="T")],
        (112086919, 112086922)
    ),
    # A stop mutation which destroys most of the protein
    (
        "9_10insTAG",
        [Variant(start=44, end=44, inserted="TAG")],
        (112086916, 112094967)
    ),
    # # A frameshift that is restored by an insertion
    # Note this gives two separate, adjacent regions
    (
        "[9_10insA;20del]",
        [
            Variant(start=44, end=44, inserted="A"),
            Variant(start=54, end=55),
        ],
        [
            (112086919, 112086925),
            (112086925, 112086928),
        ]
    ),
    # # A bigger deletion
    (
        "13_21del",
        [Variant(start=47, end=56)],
        (112086919, 112086928)
    ),
    # # An SNP that creates a STOP codon
    (
        "14G>A",
        [Variant(start=48, end=49, inserted="A", deleted="G")],
        (112086919, 112094967)
    ),
]
# fmt: on
@pytest.mark.parametrize("cdot, variants, expected", FORWARD_MUTATIONS_VARIANT)
def test_mutation_to_cds_effect_forward_new(
    cdot: str, variants: Sequence[Variant], expected: Tuple[int, int]
) -> None:
    """
    GIVEN a HGVS transcript description for a transcript on the reverse strand
    WHEN we determine the CDS effect
    THEN we should get genome coordinates
    """
    hgvs = "ENST00000375549.8:c.="
    d = Description(hgvs)
    _init_model(d)

    if not isinstance(expected, list):
        e = [expected]
    else:
        e = expected

    assert mutation_to_cds_effect(d, variants) == e


# fmt: off
MUTATIONS_VARIANT = [
    # Variant, coordinates on the genome
    # A simple missense that changes a single amino acids
    (
        "13T>A",
        [Variant(start=47573, end=47574, inserted="A", deleted="T")],
        (32435345, 32435348)
    ),
    # A stop mutation which destroys most of the protein
    (
        "9_10insTAG",
        [Variant(start=47577, end=47577, inserted="TAG")],
        (32389060, 32435351)
    ),
    # # # A frameshift that is restored by an insertion
    (
        "[10del;20_21insA]",
        [
            Variant(start=47576, end=47577),
            Variant(start=47566, end=47566, inserted="A"),
        ],
        (32435339, 32435351)
    ),
    # # # A frameshift that is restored by a bigger insertion
    (
        "[10del;20_21insATCGAATATGGGG]",
        [
            Variant(start=47566, end=47566, inserted="ATCGAATATGGGG"),
            Variant(start=47576, end=47577),
        ],
        (32435339, 32435351)),
    # # # A bigger deletion
    (
        "11_19del",
        [Variant(start=47567, end=47576)],
         (32435342, 32435351)
    ),
    # # # An inframe deletion that creates a STOP codon
    (
        "87_89del",
        [Variant(start=47497, end=47500)],
        (32389060, 32435276)
    ),
]
# fmt: on
@pytest.mark.parametrize("cdot, variants, expected", MUTATIONS_VARIANT)
def test_mutation_to_cds_effect_reverse_new(
    cdot: str, variants: Sequence[Variant], expected: Tuple[int, int]
) -> None:
    """
    GIVEN a HGVS transcript description for a transcript on the reverse strand
    WHEN we determine the CDS effect
    THEN we should get genome coordinates
    """
    hgvs = f"ENST00000452863.10:c.="
    d = Description(hgvs)
    _init_model(d)

    assert mutation_to_cds_effect(d, variants) == [expected]


RT_VARIANTS = [
    "10C>T",
    "10del",
    "10_11insA",
    "10_11delinsGG",
]


@pytest.mark.parametrize("variant", RT_VARIANTS)
def test_Variant_hgvs_round_trip_forward(
    SDHD_description: Description, variant: str
) -> None:
    """Test converting between mutalyzer delins model and Variant"""
    delins_model = _cdot_to_internal_delins(SDHD_description, CdotVariant(variant))[0]
    v = Variant.from_model(delins_model)
    assert v.to_model() == delins_model


TO_HGVS = [
    # SNP
    (Variant(44, 45, "T", "C"), "10C>T"),
    # Deletion
    (Variant(44, 45), "10del"),
    # Insertion
    (Variant(45, 45, "A"), "10_11insA"),
    # Insertion/Deletion
    (Variant(44, 46, "GG"), "10_11delinsGG"),
]


@pytest.mark.parametrize("variant, expected", TO_HGVS)
def test_Variant_to_hgvs(
    SDHD_description: Description, variant: Variant, expected: str
) -> None:
    assert to_cdot_hgvs(SDHD_description, [variant]) == expected
