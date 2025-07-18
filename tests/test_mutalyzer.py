from collections.abc import Sequence
import pytest
from pathlib import Path
from gtgt.mutalyzer import (
    CdotVariant,
    HGVS_to_genome_range,
    InternalVariant,
    append_mutation,
    combine_variants_deletion,
    exonskip,
    changed_protein_positions,
    get_exons,
    mutation_to_cds_effect,
    _cdot_to_internal_delins,
    _init_model,
    _Variant,
    to_cdot_hgvs,
)
from gtgt.mutalyzer import HGVS, Variant, variant_to_model
from mutalyzer.description import Description
from mutalyzer.converter.to_delins import variants_to_delins
from mutalyzer.converter.to_internal_coordinates import to_internal_coordinates
from mutalyzer.converter.to_internal_indexing import to_internal_indexing
from mutalyzer.description_model import get_reference_id, variants_to_description
from gtgt.transcript import Transcript
from gtgt.models import TranscriptModel
import json
import copy

from itertools import zip_longest
from typing import Any, List, Tuple, Dict


def retrieve_raw(
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


def get_cds_to_mrna(cds_id: Any, timeout: Any = 10) -> None:
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
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", retrieve_raw)
    monkeypatch.setattr("mutalyzer.description.get_cds_to_mrna", get_cds_to_mrna)


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


POSITIONS = [
    # c., genome range
    # Start of the transcript / first exon
    ("-35", (0, 1)),
    # Start of the coding region
    ("1", (35, 36)),
    # End of the first exon
    ("52", (86, 87)),
    # Start of the second exon
    ("53", (984, 985)),
    # End of the second exon
    ("169", (1100, 1101)),
    # Start of the third exon
    ("170", (1994, 1995)),
    # End of the third exon
    ("314", (2138, 2139)),
    # Start of the fourth exon
    ("315", (7932, 7933)),
    # End of the CDS
    ("480", (8097, 8098)),
    # End of the transcript / fourth exon
    ("*824", (8921, 8922)),
]


@pytest.mark.parametrize("cdot, genomic", POSITIONS)
def test_genomic_to_c_dot_SDHD(cdot: str, genomic: Tuple[int, int]) -> None:
    """
    GIVEN a Transcript and a position in c. notation
    WHEN we translate the c. to genomic positions
    THEN we should get a range, on the genomic coordinates

    NOTE that genomic coordinates refer to the UCSC annotations on the genome,
    i.e. 0 based, half open. Not to be confused with hgvs g. positions
    """
    # Ensemble transcript ID for SDHD
    ENST = "ENST00000375549.8"
    # Build the HGVS description
    description = f"{ENST}:c.{cdot}"

    d = Description(description)
    _init_model(d)

    start, end = HGVS_to_genome_range(d)

    # Genomic offset of SDHD, just used for testing so the genomic positions
    # are manageable
    offset = 112086872
    start -= offset
    end -= offset

    assert (start, end) == genomic


POSITIONS = [
    # c., genome range
    # Start of the transcript / first exon
    ("-179", (47764, 47765)),
    # Start of the coding region
    ("1", (47585, 47586)),
    # End of the first exon
    ("661", (46925, 46926)),
    # Start of the second exon
    ("662", (40844, 40845)),
    # End of the second exon
    ("784", (40722, 40723)),
    # Start of the third exon
    ("785", (40283, 40284)),
    # End of the third exon
    ("887", (40181, 40182)),
    # Start of the tenth exon
    ("1448", (1404, 1405)),
    # End of the CDS
    ("1569", (1283, 1284)),
    # End of the transcript / tenth exon
    ("*1283", (0, 1)),
]


@pytest.mark.parametrize("cdot, genomic", POSITIONS)
def test_genomic_to_c_dot_WT1(cdot: str, genomic: Tuple[int, int]) -> None:
    """
    GIVEN a Transcript and a position in c. notation
    WHEN we translate the c. to genomic positions
    THEN we should get a range, on the genomic coordinates

    NOTE: WT1 is on the reverse strand

    NOTE that genomic coordinates refer to the UCSC annotations on the genome,
    i.e. 0 based, half open. Not to be confused with hgvs g. positions
    """
    # Ensemble transcript ID for SDHD
    ENST = "ENST00000452863.10"
    # Build the HGVS description
    description = f"{ENST}:c.{cdot}"

    d = Description(description)
    _init_model(d)

    start, end = HGVS_to_genome_range(d)

    # Genomic offset of WT1, just used for testing so the genomic positions are
    # manageable
    offset = 32387774
    start -= offset
    end -= offset

    assert (start, end) == genomic


SUPPORTED_DESCRIPTIONS = [
    # SNP
    ("ENST00000375549.8:c.5C>T"),
    # Deletion of range
    ("ENST00000375549.8:c.5_15del"),
    # Inversion
    ("ENST00000375549.8:c.5_15inv"),
    # Empty position
    ("ENST00000375549.8:c.5"),
    # Empty range
    ("ENST00000375549.8:c.5_10"),
    # Insertion of a range
    ("ENST00000375549.8:c.10_11ins50_60"),
]


@pytest.mark.parametrize("description", SUPPORTED_DESCRIPTIONS)
def test_supported_descriptions(description: str) -> None:
    d = Description(description)
    _init_model(d)
    HGVS_to_genome_range(d)


UNSUPPORTED_DESCRIPTIONS = [
    # Description with multiple variants
    ("ENST00000375549.8:c.[5C>T;10del]"),
    # Empty HGVS description
    ("ENST00000375549.8:c.="),
    # Delins of a range
    ("ENST00000375549.8:c.10_15delins50_60"),
]


@pytest.mark.parametrize("description", UNSUPPORTED_DESCRIPTIONS)
def test_unsupported_descriptions(description: str) -> None:
    with pytest.raises(ValueError):
        d = Description(description)
        _init_model(d)
        HGVS_to_genome_range(d)


def test_exonskip_SDHD() -> None:
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    results = [
        "ENST00000375549.8:c.53_169del",
        "ENST00000375549.8:c.170_314del",
    ]
    for output, expected in zip_longest(exonskip(d), results):
        assert output.hgvs == expected


def test_new_exonskip_SDHD() -> None:
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    results = [
        "ENST00000375549.8:c.53_169del",
        "ENST00000375549.8:c.170_314del",
    ]
    for output, expected in zip_longest(exonskip(d), results):
        assert output.hgvs == expected


def test_exonskip_WT1() -> None:
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
    for output, expected in zip_longest(exonskip(d), results):
        assert output.hgvs == expected


MUTATIONS = [
    # HGVS, coordinates on the genome,
    # A simple missense that changes a single amino acids
    ("13T>A", (32435345, 32435348)),
    # A stop mutation which destroys most of the protein
    ("9_10insTAG", (32389060, 32435351)),
    # # A frameshift that is restored by an insertion
    ("[10del;20_21insA]", (32435339, 32435351)),
    # # A frameshift that is restored by a bigger insertion
    ("[10del;20_21insATCGAATATGGGG]", (32435339, 32435351)),
    # # A bigger deletion
    ("11_19del", (32435342, 32435351)),
    # # An inframe deletion that creates a STOP codon
    ("87_89del", (32389060, 32435276)),
]


@pytest.mark.parametrize("variants, expected", MUTATIONS)
def test_mutation_to_cds_effect(
    variants: CdotVariant, expected: Tuple[int, int]
) -> None:
    """
    GIVEN a HGVS transcript description for a transcript on the reverse strand
    WHEN we determine the CDS effect
    THEN we should get genome coordinates
    """
    d = Description("ENST00000452863.10:c.=")
    _init_model(d)

    assert mutation_to_cds_effect(d, variants) == [expected]


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


VARIANTS = [
    # variant, Transcript effect
    # In frame deletion that creates a STOP codon
    ("ENST00000452863.10:c.87_89del", 0.0018),
    # In frame deletion that does not make a STOP
    ("ENST00000452863.10:c.85_87del", 0.9999),
    # Synonymous mutation
    ("ENST00000452863.10:c.13T>C", 1),
]


@pytest.mark.parametrize("hgvs_description, effect", VARIANTS)
def test_mutate_transcript_with_variant(
    hgvs_description: str, effect: float, WT: Transcript
) -> None:

    transcript_id = hgvs_description.split(":c.")[0]
    variants = CdotVariant(hgvs_description.split(":c.")[1])

    d = Description(f"{transcript_id}:c.=")
    _init_model(d)

    modified = copy.deepcopy(WT)
    modified.mutate(d, variants)

    cmp = modified.compare(WT)

    coding_exons = cmp[1]
    assert coding_exons.name == "coding_exons"
    assert coding_exons.percentage == pytest.approx(effect, abs=0.0001)


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
def test_variant_to_model(variant: CdotVariant, variant_models: List[Variant]) -> None:
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
    v = _Variant(10, 11, "ATG")
    assert str(v) == "Variant(start=10, end=11, inserted=ATG, deleted=)"


def test_Variant_class_str_snp() -> None:
    """SNPS are special, since they contain the inserted sequence"""
    # 10A>T
    v = _Variant(10, 11, "T", "A")
    assert str(v) == "Variant(start=10, end=11, inserted=T, deleted=A)"


def test_Variant_class_to_model_positions() -> None:
    """Test converting a variant to model"""
    v = _Variant(10, 11, "ATG")
    model = v.to_model()

    assert model["location"]["start"]["position"] == 10
    assert model["location"]["end"]["position"] == 11
    assert model["inserted"][0]["sequence"] == "ATG"

    # Deleted entry is missing for deletions
    assert "deleted" not in model


def test_Variant_class_to_model_snp() -> None:
    # 10 A>T
    v = _Variant(10, 11, "T", "A")
    model = v.to_model()

    assert model["inserted"][0]["sequence"] == "T"
    assert model["deleted"][0]["sequence"] == "A"


def test_Variant_deleted_snp_only() -> None:
    """Test that deleted is only defined for SNPs, not larger indels

    (To match the behaviour of Mutalyzer)
    """
    with pytest.raises(ValueError):
        # 10_12delinsGG
        _Variant(10, 12, "GG", "AT")


def test_Variant_class_no_inserted_sequence() -> None:
    """
    GIVEN a Variant with an empty inserted sequence
    WHEN we convert to a model
    THEN inserted should be an empty list
    """
    v = _Variant(10, 11, "")
    model = v.to_model()

    assert model["inserted"] == []


def test_Variant_class_end_after_start() -> None:
    """Ensure that end is after start"""
    with pytest.raises(ValueError):
        _Variant(11, 10, "")


ORDERING = [
    # Ends are touching
    (_Variant(1, 3), _Variant(3, 5)),
    # Gap between variants
    (_Variant(0, 1), _Variant(2, 4)),
]


@pytest.mark.parametrize("a, b", ORDERING)
def test_Variant_class_relative_positions(a: _Variant, b: _Variant) -> None:
    """Variant a is before variant b"""
    assert a.before(b)
    assert b.after(a)


INSIDE = [
    # Variants are inside themselves
    (_Variant(0, 3), _Variant(0, 3)),
    # Smaller variant a is inside b
    (_Variant(1, 3), _Variant(0, 3)),
    (_Variant(0, 2), _Variant(0, 3)),
    (_Variant(1, 2), _Variant(0, 3)),
]


@pytest.mark.parametrize("a, b", INSIDE)
def test_Variant_class_inside(a: _Variant, b: _Variant) -> None:
    """Variant a is inside variant b"""
    assert a.inside(b)


NOT_INSIDE = [
    # a starts outside of b
    (_Variant(0, 3), _Variant(1, 3)),
    # a ends outside of b
    (_Variant(1, 4), _Variant(1, 3)),
    # b is inside of a
    (_Variant(0, 3), _Variant(1, 2)),
    # a before b
    (_Variant(0, 3), _Variant(3, 5)),
    # a after b
    (_Variant(3, 5), _Variant(1, 3)),
]


@pytest.mark.parametrize("a, b", NOT_INSIDE)
def test_Variant_class_not_inside(a: _Variant, b: _Variant) -> None:
    """Variant a is not inside variant b"""
    assert not a.inside(b)


OVERLAP = [
    # b ends inside a
    (_Variant(2, 5), _Variant(1, 3)),
    # b fully inside a
    (_Variant(2, 5), _Variant(3, 4)),
    # b fully inside a, ends at end
    (_Variant(2, 5), _Variant(3, 5)),
    # b ends inside a
    (_Variant(2, 5), _Variant(2, 4)),
    # b starts in a, extends after
    (_Variant(2, 5), _Variant(3, 6)),
    # b start before a
    (_Variant(2, 5), _Variant(1, 5)),
    # a is inside b
    (_Variant(2, 5), _Variant(1, 6)),
    # a is inside b
    (_Variant(2, 5), _Variant(2, 6)),
    # a equals b
    (_Variant(2, 5), _Variant(2, 5)),
]


@pytest.mark.parametrize("a, b", OVERLAP)
def test_Variant_class_overlap(a: _Variant, b: _Variant) -> None:
    """Variant a and b overlap"""
    assert a.overlap(b)
    assert b.overlap(a)


NO_OVERLAP = [
    # a is after b
    (_Variant(2, 5), _Variant(1, 2)),
    # a is before b
    (_Variant(2, 5), _Variant(5, 6)),
]


@pytest.mark.parametrize("a, b", NO_OVERLAP)
def test_Variant_class_no_overlap(a: _Variant, b: _Variant) -> None:
    """Variant a and b do not overlap"""
    assert not a.overlap(b)
    assert not b.overlap(a)


def test_Variant_class_order() -> None:
    """Test sorting variants by start position"""
    v1 = _Variant(10, 11)
    v2 = _Variant(0, 1)
    assert sorted([v1, v2]) == [v2, v1]


def test_Variant_class_order_error() -> None:
    """Test error when sorting overlapping variants"""
    v1 = _Variant(10, 11)
    v2 = _Variant(10, 11)
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
    v = _Variant.from_model(delins_model)
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
    v = _Variant.from_model(delins_model)
    assert v.start == 0
    assert v.end == 2
    # Test that the emtpy sequence is a string, not a list
    assert v.inserted == ""


def test_combine_variants_deletion_empty() -> None:
    variants: List[_Variant] = list()
    deletion = _Variant(0, 10)
    assert combine_variants_deletion(variants, deletion) == [deletion]


# fmt: off
COMBINE = [
    ( # The variants are out of order
        # Variants
        [_Variant(5, 7), _Variant(2, 4)],
        # Deletion
        _Variant(10, 11),
        # Expected
        [_Variant(2,4), _Variant(5, 7), _Variant(10, 11)],
    ),
    ( # The deletion is before both variants
        # Variants
        [_Variant(2, 4), _Variant(5, 7)],
        # Deletion
        _Variant(0, 1),
        # Expected
        [_Variant(0, 1), _Variant(2,4), _Variant(5, 7)],
    ),
    ( # The deletion is between the variants
        # Variants
        [_Variant(2, 4), _Variant(5, 7)],
        # Deletion
        _Variant(4, 5),
        # Expected
        [_Variant(2,4), _Variant(4, 5), _Variant(5, 7)],
    ),
    ( # The deletion is after both variants
        # Variants
        [_Variant(2, 4), _Variant(5, 7)],
        # Deletion
        _Variant(10, 11),
        # Expected
        [_Variant(2,4), _Variant(5, 7), _Variant(10, 11)],
    ),
    ( # The deletion contains the first variant
        # Variants
        [_Variant(2, 4), _Variant(5, 7)],
        # Deletion
        _Variant(1, 5),
        # Expected
        [_Variant(1, 5), _Variant(5, 7)],
    ),
    ( # The deletion contains the second variant
        # Variants
        [_Variant(2,4), _Variant(5, 7)],
        # Deletion
        _Variant(4, 11),
        # Expected
        [_Variant(2,4), _Variant(4, 11)],
    ),
    ( # The deletion contains both variants
        # Variants
        [_Variant(2,4), _Variant(5, 7)],
        # Deletion
        _Variant(2, 7),
        # Expected
        [_Variant(2, 7)],
    ),
]
# fmt: on


@pytest.mark.parametrize("variants, deletion, expected", COMBINE)
def test_combine_variants_deletion(
    variants: Sequence[_Variant], deletion: _Variant, expected: Sequence[_Variant]
) -> None:
    combined = combine_variants_deletion(variants, deletion)
    assert combined == expected


def test_combine_variants_deletion_variants_overlap_eachother() -> None:
    """Test that we get a value error if the variants overlap"""
    variants = [_Variant(0, 2), _Variant(1, 3)]
    deletion = _Variant(10, 11)
    with pytest.raises(ValueError):
        combine_variants_deletion(variants, deletion)


def test_combine_variants_deletion_variants_partially_overlap_deletion() -> None:
    """Test that we get a value error if one the variants partially overlaps the deletion"""
    variants = [_Variant(2, 4)]
    deletion = _Variant(3, 11)
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
    v = _Variant.from_model(delins_model)
    assert v.to_model() == delins_model


TO_HGVS = [
    # SNP
    (_Variant(44, 45, "T", "C"), "10C>T"),
    # Deletion
    (_Variant(44, 45), "10del"),
    # Insertion
    (_Variant(44, 45, "A"), "10_11insA"),
    # Insertion/Deletion
    (_Variant(44, 46, "GG"), "10_11delinsGG"),
]


# @pytest.mark.parametrize("variant, expected", TO_HGVS)
# def test_Variant_to_hgvs(
#     SDHD_description: Description, variant: _Variant, expected: str
# ) -> None:
#     print()
#     assert to_cdot_hgvs(SDHD_description, [variant]) == expected
