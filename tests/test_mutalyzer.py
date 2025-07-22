from collections.abc import Sequence
import pytest
from pathlib import Path
from gtgt.mutalyzer import (
    CdotVariant,
    append_mutation,
    cds_to_internal_positions,
    combine_variants_deletion,
    changed_protein_positions,
    get_exons,
    mutation_to_cds_effect,
    _cdot_to_internal_delins,
    _init_model,
    Variant,
    skip_adjacent_exons,
    to_cdot_hgvs,
    sliding_window,
    _exon_string,
)
from gtgt.mutalyzer import Variant_Dict, variant_to_model
from mutalyzer.description import Description
from mutalyzer.description_model import variants_to_description
from gtgt.transcript import Transcript
from gtgt.models import TranscriptModel
import json

from itertools import zip_longest
from typing import Any, List, Tuple


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


def test_one_adjacent_exonskip_SDHD() -> None:
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    results = [
        "ENST00000375549.8:c.53_169del",
        "ENST00000375549.8:c.170_314del",
    ]
    for output, expected in zip_longest(skip_adjacent_exons(d), results):
        assert output.hgvs == expected


def test_two_adjacent_exonskip_SDHD() -> None:
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    results = [
        "ENST00000375549.8:c.53_314del",
    ]
    for output, expected in zip_longest(
        skip_adjacent_exons(d, number_to_skip=2), results
    ):
        assert output.hgvs == expected


def test_no_possible_exonskip_SDHD() -> None:
    """
    GIVEN a transcript with 4 exons (2 can be skipped)
    WHEN we try to skip 3 adjacent exons
    THEN we should get an empty list of therapies
    """
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    assert skip_adjacent_exons(d, number_to_skip=3) == list()


def test_one_adjacent_exonskip_WT1() -> None:
    d = Description("ENST00000452863.10:c.=")
    _init_model(d)
    results = [
        "ENST00000452863.10:c.662_784del",
        "ENST00000452863.10:c.785_887del",
        "ENST00000452863.10:c.888_965del",
        "ENST00000452863.10:c.966_1016del",
        "ENST00000452863.10:c.1017_1113del",
        "ENST00000452863.10:c.1114_1264del",
        "ENST00000452863.10:c.1265_1354del",
        "ENST00000452863.10:c.1355_1447del",
    ]
    # for output, expected in zip_longest(skip_adjacent_exons(d), results):
    for output, expected in zip_longest(skip_adjacent_exons(d), results):
        assert output.hgvs == expected


def test_two_adjacent_exonskip_WT1() -> None:
    d = Description("ENST00000452863.10:c.=")
    _init_model(d)
    results = [
        "ENST00000452863.10:c.662_887del",
        "ENST00000452863.10:c.785_965del",
        "ENST00000452863.10:c.888_1016del",
        "ENST00000452863.10:c.966_1113del",
        "ENST00000452863.10:c.1017_1264del",
        "ENST00000452863.10:c.1114_1354del",
        "ENST00000452863.10:c.1265_1447del",
    ]
    for output, expected in zip_longest(skip_adjacent_exons(d, 2), results):
        assert output.hgvs == expected


def test_sliding_window_size_one() -> None:
    s = "ABCDEF"
    # assert list(sliding_window(s, 1)) == [[x] for x in "A B C D E F".split()]
    assert list(sliding_window(s, 1)) == [["A"], ["B"], ["C"], ["D"], ["E"], ["F"]]


def test_sliding_window_size_2() -> None:
    s = "ABCDEF"
    assert list(sliding_window(s, 2)) == [
        ["A", "B"],
        ["B", "C"],
        ["C", "D"],
        ["D", "E"],
        ["E", "F"],
    ]


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


CDOT_MUTATIONS = [
    # c. variant, i. delins variant
    ("13T>A", "47_48delTinsA"),
    ("-35_52del", "0_87delins"),
    ("53_169del", "984_1101delins"),
    ("315_*824del", "7932_8922delins"),
    # Multiple variants
    ("[13T>A;15_16delinsAT]", "[47_48delTinsA;49_51delinsAT]"),
]


@pytest.mark.parametrize("cdot, internal_delins", CDOT_MUTATIONS)
def test_cdot_to_indel(cdot: CdotVariant, internal_delins: str) -> None:
    """
    GIVEN a list of c. mutations a string
    WHEN we convert them to internal i. indels
    THEN we should get a list of internal variants
    """
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    indels = _cdot_to_internal_delins(d, cdot)
    assert variants_to_description(indels) == internal_delins


@pytest.fixture
def WT() -> Transcript:
    """
    Transcript for WT1, using real genomic positions
    """
    path = "tests/data/ENST00000452863.10.Transcript.json"
    with open(path) as fin:
        js = json.load(fin)

    t = TranscriptModel.model_validate(js)

    return t.to_transcript()


def test_analyze_transcript(WT: Transcript) -> None:
    # In frame deletion that creates a STOP codon
    # variant = "ENST00000452863.10:c.87_89del"
    # Frameshift in small in-frame exon 5
    variant = "ENST00000452863.10:c.970del"

    results = WT.analyze(variant)

    wildtype = results[0]
    coding_exons = wildtype.comparison[1]
    assert coding_exons.name == "coding_exons"
    assert coding_exons.percentage == 1.0


def test_analyze_transcript_r_coordinate(WT: Transcript) -> None:
    """Test analyzing a transcript using the r. coordinate system"""
    # In frame deletion that creates a STOP codon
    variant = "ENST00000452863.10:r.970del"

    results = WT.analyze(variant)

    wildtype = results[0]
    coding_exons = wildtype.comparison[1]
    assert coding_exons.name == "coding_exons"
    assert coding_exons.percentage == 1.0


APPEND_VARIANT = [
    # mutation, predicted protein description as readout
    ("10del", "MAVSGG*")
]


@pytest.mark.parametrize("mutation, protein", APPEND_VARIANT)
def test_append_mutation(mutation: CdotVariant, protein: str) -> None:
    """
    GIVEN a string denoting a HGVS variant
    WHEN we append this variant to an existing Description
    THEN the Description should be updated
    """
    transcript = "ENST00000375549.8:c.="
    d = Description(f"{transcript}")
    _init_model(d)

    append_mutation(d, mutation)

    assert d.output()["protein"]["predicted"] == protein


APPEND_TO_EXISTING = [
    # Add an insertion after the deletion to restore the reading frame
    ("10del", "11_12insA", "MAVYWRLSAV"),
    # Add an insertion before the deletion, e.g. the variants are out of order after appending
    ("10del", "4dup", "MGGFWRLSAV"),
    # Existing variant is a snp, not deletion
    ("10C>T", "11dup", "MAVFLEAECR"),
    # Existing variant is a dup
    ("10dup", "15del", "MAVPLRLSAV"),
    # Existing variant is an inversion
    ("10_15inv", "16A>T", "MAVPEWLSAV"),
    # Existing variant is an delins
    ("10_12delinsGG", "15dup", "MAVGGRLSAV"),
    # 8del gets normalized to 9del
    ("8del", "9del", "MAALEAECRL"),
    ("9del", "8del", "MAALEAECRL"),
]


@pytest.mark.parametrize("existing, novel, protein", APPEND_TO_EXISTING)
def test_append_mutation_to_existing_variant(
    existing: CdotVariant, novel: CdotVariant, protein: str
) -> None:
    """
    GIVEN a string denoting a HGVS variant
    WHEN we append this variant to an existing Description
    THEN the Description should be updated
    """
    transcript = f"ENST00000375549.8:c.{existing}"
    d = Description(transcript)
    _init_model(d)
    append_mutation(d, novel)

    # Test the first 10 aa of the protein as readout
    assert d.output()["protein"]["predicted"][:10] == protein


APPEND_OVERLAPPING_VARIANT = [
    # Transcript variant, new variant
    # Add another variant at the same location
    ("10C>T", "10del"),
    # Both variants we want to add overlap
    ("=", "[10del;10C>T]"),
    # Both variants already in the description overlap
    ("[10del;10C>T]", "15del"),
    # Overlap with an inversion
    ("10_15inv", "12C>T"),
    ("10_15inv", "12G>T"),
    ("10_15inv", "10C>T"),
]


@pytest.mark.parametrize("existing, novel", APPEND_OVERLAPPING_VARIANT)
def test_appending_overlapping_variants(
    existing: CdotVariant, novel: CdotVariant
) -> None:
    """
    GIVEN a transcript with existing variant(s)
    WHEN we attempt to add novel variant which overlaps the existing variants
    THEN we should throw a ValueError
    """
    transcript = f"ENST00000375549.8:c.{existing}"
    d = Description(transcript)
    _init_model(d)

    with pytest.raises(ValueError):
        append_mutation(d, novel)


PARSE_VARIANT = [
    (
        "10del",
        [
            {
                "location": {"type": "point", "position": 10},
                "type": "deletion",
                "source": "reference",
            }
        ],
    ),
    (
        "[10del; 16dup]",
        [
            {
                "location": {"type": "point", "position": 10},
                "type": "deletion",
                "source": "reference",
            },
            {
                "location": {"type": "point", "position": 16},
                "type": "duplication",
                "source": "reference",
            },
        ],
    ),
]


@pytest.mark.parametrize("variant, variant_models", PARSE_VARIANT)
def test_variant_to_model(
    variant: CdotVariant, variant_models: List[Variant_Dict]
) -> None:
    """
    GIVEN a string denoting a HGVS variant
    WHEN we parse this into a Variant
    THEN it should contain the expected values
    """
    assert variant_to_model(variant) == variant_models


PROTEIN_EXTRACTOR = [
    # No protein description
    ("", "", []),
    # No change
    ("A", "A", []),
    # A single change
    ("A", "T", [(0, 1)]),
    # A single change on the second position
    ("AAA", "ATA", [(1, 2)]),
    # Change in a repeat region
    ("AA", "A", [(1, 2)]),
    ("AAA", "A", [(1, 3)]),
    # A delins
    ("AAA", "ATC", [(1, 3)]),
    # An insertion, which we ignore
    ("AAA", "AATA", []),
    # A delins of TAG, which is equivalent to two insertions
    ("AAA", "ATAGAA", []),
    # A delins which is equivalent to a deletion
    ("AAA", "AGGGA", [(1, 2)]),
    # Multiple deletions
    ("AAA", "TAT", [(0, 1), (2, 3)]),
]


@pytest.mark.parametrize("reference, observed, expected", PROTEIN_EXTRACTOR)
def test_changed_protein_positions(
    reference: str, observed: str, expected: List[Tuple[int, int]]
) -> None:
    """
    GIVEN a referene and observed sequence
    WHEN we extrat the protein changes
    THEN we should get 0 based positions of the deleted residues
    """
    assert changed_protein_positions(reference, observed) == expected


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


def test_get_exons_forward(SDHD_description: Description) -> None:
    """Text extracting exons from a Description object"""
    expected = (0, 87)

    assert get_exons(SDHD_description, in_transcript_order=True)[0] == expected
    assert get_exons(SDHD_description, in_transcript_order=False)[0] == expected


def test_exons_forward(WT1_description: Description) -> None:
    """Text extracting exons from a Description object"""
    expected_transcript_order = (46925, 47765)
    expected_genomic_order = (0, 1405)

    assert (
        get_exons(WT1_description, in_transcript_order=True)[0]
        == expected_transcript_order
    )
    assert (
        get_exons(WT1_description, in_transcript_order=False)[0]
        == expected_genomic_order
    )


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


# fmt: off
CDS_POSITIONS = [
    (0, 0),
    (4, 4),
    (5, 11),
    (7, 13),
    (8, 100)
]
# fmt: on


@pytest.mark.parametrize("position, expected", CDS_POSITIONS)
def test_cds_to_internal_position_forward(position: int, expected: int) -> None:
    exons = [(0, 5), (11, 14), (100, 120)]
    assert cds_to_internal_positions(position, exons) == expected


CDS_POSITIONS_REV = [
    (0, 13),
    (2, 11),
    (3, 4),
    (7, 0),
]


@pytest.mark.parametrize("position, expected", CDS_POSITIONS_REV)
def test_cds_to_internal_position_reverse(position: int, expected: int) -> None:
    exons = [(0, 5), (11, 14)]
    assert cds_to_internal_positions(position, exons, reverse=True) == expected


CDS_POSITIONS_OFFSET = [
    # internal, forward, reverse
    (0, 12, 39),
    (1, 25, 38),
    (2, 26, 37),
    (3, 27, 29),
    (5, 29, 27),
    (6, 37, 26),
    (7, 38, 25),
    (8, 39, 12),
]


@pytest.mark.parametrize("position, forward, reverse", CDS_POSITIONS_OFFSET)
def test_cds_to_internal_position(position: int, forward: int, reverse: int) -> None:
    # coding exons that do not start at position 0
    # Size:      1        4         2
    exons = [(12, 13), (25, 30), (37, 40)]
    assert cds_to_internal_positions(position, exons, reverse=False) == forward
    assert cds_to_internal_positions(position, exons, reverse=True) == reverse


def test_cds_to_internal_positions_out_of_range() -> None:
    exons = [(0, 5)]

    # position 4 is still in the exon
    cds_to_internal_positions(4, exons)
    cds_to_internal_positions(4, exons, reverse=True)

    # position 5 is outside the exon
    with pytest.raises(ValueError):
        cds_to_internal_positions(5, exons)
    with pytest.raises(ValueError):
        cds_to_internal_positions(5, exons, reverse=True)

    # CDS that does not start at o
    exons = [(12, 13), (25, 30), (37, 40)]

    with pytest.raises(ValueError):
        cds_to_internal_positions(9, exons)
    with pytest.raises(ValueError):
        cds_to_internal_positions(9, exons, reverse=True)


EXON_DESCRIPTION = [
    ([2], "exon 2"),
    ([3, 5], "exons 3 and 5"),
    ([3, 4, 5, 6], "exons 3, 4, 5 and 6"),
]


@pytest.mark.parametrize("exons, expected", EXON_DESCRIPTION)
def test_exon_string(exons: Sequence[int], expected: str) -> None:
    assert _exon_string(exons) == expected
