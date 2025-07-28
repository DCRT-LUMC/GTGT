#!/usr/bin/env python3

import json
import sys

import pytest
from mutalyzer.description import Description
from mutalyzer.description_model import get_reference_id, variants_to_description
from mutalyzer.protein import get_protein_description
from mutalyzer.reference import get_protein_selector_model

from gtgt.mutalyzer import Variant


def pprint(thing):
    print(json.dumps(thing, indent=True, sort_keys=True))


def _init_d(hgvs_description):
    """
    No normalization is performed! The extractor is 'short-circuited'.
    """
    d = Description(hgvs_description)
    d.to_delins()
    d.de_hgvs_internal_indexing_model = d.delins_model
    d.construct_de_hgvs_internal_indexing_model()
    d.construct_de_hgvs_coordinates_model()
    d.construct_normalized_description()
    d.construct_protein_description()
    return d


def delins_from_description(d):
    """Return the internal delins model from a Description object"""
    delins_models = d.de_hgvs_internal_indexing_model["variants"]
    assert len(delins_models) == 1
    return d.de_hgvs_internal_indexing_model["variants"][0]


def protein_from_description(d):
    """Return the predicted protein sequence form a description"""
    return d.protein["predicted"]


def sequence_from_description(d):
    """Return the sequence form a description"""
    _id = d.input_model["reference"]["id"]
    return d.references[_id]["sequence"]["seq"]


def protein_from_variant(v, d):
    """Return the predicted protein sequence from a Variant"""

    ref_id = get_reference_id(d.corrected_model)

    selector_model = get_protein_selector_model(
        d.references[ref_id]["annotations"], ref_id
    )

    # Get the sequence (needed for complex variants)
    delins_models = [v.to_model()]

    desc, ref, predicted, *rest = get_protein_description(
        delins_models, d.references, selector_model
    )
    return predicted


class TestForward:
    """Test equivalent models on the forward strand"""

    # SDHD, forward strand
    transcript = "ENST00000375549.8"
    empty_transcript = _init_d(f"{transcript}:c.=")

    VARIANTS = [
        # A SNP
        ("10C>T", Variant(44, 45, inserted="T", deleted="C")),
        # A deletion
        ("10del", Variant(44, 45, inserted="")),
        # An insertion
        ("10_11insA", Variant(45, 45, inserted="A")),
        # Delins version of 10C>T (Note that the deleted part is lost)
        ("10_10delinsT", Variant(44, 45, inserted="T")),
        # A duplication
        ("10dup", Variant(45, 45, inserted="C")),
        # The same duplication, with Variant as a delins (the deleted part is
        # implicit)
        ("10dup", Variant(44, 45, inserted="CC")),
        # A duplication
        ("10_11dup", Variant(44, 44, inserted="CT")),
        # The same duplication, where Variant deletes the first "CT", and then
        # inserts it twice
        ("10_11dup", Variant(44, 46, inserted="CTCT")),
        # Inversion, equivalent to 10C>G
        ("10_10inv", Variant(44, 45, inserted="G")),
        ("10C>G", Variant(44, 45, inserted="G")),
        # Inversion, equivalent to 10_11delinsAG
        ("10_11inv", Variant(44, 46, inserted="AG")),
        ("10_11delinsAG", Variant(44, 46, inserted="AG")),
        # Inversion, not symetrical
        ("18_20inv", Variant(52, 55, inserted="AGC")),
        ("18_20delinsAGC", Variant(52, 55, inserted="AGC")),
        # Small mononucleotide repeat:
        # ref: GGTT  CTCT
        # obs: GGTTTTCTCT
        # hgvs notation: 9_10insTT, 8_9T[4]
        ("8_9T[4]", Variant(44, 44, inserted="TT")),
        # Repeat of multiple distinct nucleotides
        # ref: GGTTCTCT    GGA
        # obs: GGTTCTCTCTCTGGA
        # hgvs notation: 9_10insCTCT, 10_13CT[4]
        ("10_13CT[4]", Variant(44, 44, inserted="CTCT")),
    ]

    @pytest.mark.parametrize("hgvs,variant", VARIANTS)
    def test_hgvs_Variant_equivalence_via_protein(self, hgvs: str, variant: Variant):
        """Test hgvs and Variant equivalence by comparing the protein prediction

        The goal here is to verify that the model representation of the Variant
        is usable by mutalyzer in the same way as the original HGVS description
        """
        d = _init_d(f"{self.transcript}:c.{hgvs}")
        variant_protein = protein_from_variant(variant, self.empty_transcript)
        description_protein = protein_from_description(d)

        assert variant_protein == description_protein

    VARIANTS = [
        ("10C>T", Variant(44, 45, inserted="T", deleted="C")),
        ("10del", Variant(44, 45, inserted="")),
        ("10_11insA", Variant(45, 45, inserted="A")),
        ("10_10delinsT", Variant(44, 45, inserted="T")),
        ("10dup", Variant(44, 45, inserted="CC")),
        ("10_11dup", Variant(44, 46, inserted="CTCT")),
        ("10_10inv", Variant(44, 45, inserted="G")),
        ("10_11inv", Variant(44, 46, inserted="AG")),
        ("18_20inv", Variant(52, 55, inserted="AGC")),
        # Repeat is a delins of the existing repeat units, and insertion of the
        # new repeat
        ("8_9T[4]", Variant(42, 44, inserted="TTTT")),
        ("10_13CT[4]", Variant(44, 48, inserted="CTCTCTCT")),
    ]

    @pytest.mark.parametrize("hgvs, variant", VARIANTS)
    def test_hgvs_Variant_delins_model(self, hgvs: str, variant: Variant):
        """Test hgvs and Variant delins model equivalence directly

        The goal here is to verify that the Variant.from_model structure is
        always a simple delins model, even for duplications, repeats and
        inversions.
        """
        d = _init_d(f"{self.transcript}:c.{hgvs}")
        delins_model = delins_from_description(d)

        seq = sequence_from_description(d)
        variant_model = Variant.from_model(delins_model, sequence=seq)
        assert variant_model == variant

        # assert delins_model == variant_model


class TestReverse:
    """Test equivalent models on the reverse strand"""

    # SDHD, forward strand
    transcript = "ENST00000452863.10"
    empty_transcript = _init_d(f"{transcript}:c.=")

    # Note that on the Variant the nucleotides are (reverse?) complemented
    VARIANTS = [
        # A SNP,
        ("10C>T", Variant(47576, 47577, inserted="A", deleted="G")),
        # A deletion
        ("10del", Variant(47576, 47577, inserted="")),
        # An insertion
        ("10_11insA", Variant(47576, 47576, inserted="T")),
        # Delins version of 10C>T (Note that the deleted part is lost)
        ("10_10delinsT", Variant(47576, 47577, inserted="A")),
        # A duplication
        ("10dup", Variant(47576, 47576, inserted="G")),
        # The same duplication, with Variant as a delins (the deleted part is
        # implicit)
        ("10dup", Variant(47575, 47576, inserted="GG")),
        ("10_10delinsCC", Variant(47576, 47577, inserted="GG")),
        # A duplication, note that "AG" is the reverse complement of "CT"
        ("10_11dup", Variant(47575, 47575, inserted="AG")),
        ("11_12insCT", Variant(47575, 47575, inserted="AG")),
        ("12_13delinsCTCT", Variant(47575, 47577, inserted="AGAG")),
        ("12_13delinsCTCT", Variant(47573, 47575, inserted="AGAG")),
        # Inversion, equivalent to 10C>G
        ("10_10inv", Variant(47576, 47577, inserted="C")),
        ("10C>G", Variant(47576, 47577, inserted="C")),
        # Inversion, equivalent to 10_11delinsAG
        ("10_11inv", Variant(47575, 47577, inserted="CT")),
        ("10_11delinsAG", Variant(47575, 47577, inserted="CT")),
        # Inversion, not symetrical
        ("18_20inv", Variant(47566, 47569, inserted="GCA")),
        ("18_20delinsTGC", Variant(47566, 47569, inserted="GCA")),
        # Small mononucleotide repeat
        # 7_8T[4],  Variant representation = 7_8delinsTTTT
        ("7_8T[4]", Variant(47578, 47580, inserted="AAAA")),
        # Repeat of multiple distinct nucleotides
        # 10_13CT[4], Variant representation = 10_13delinsCTCTCTCT
        ("10_13CT[4]", Variant(47573, 47577, inserted="AGAGAGAG")),
    ]

    @pytest.mark.parametrize("hgvs,variant", VARIANTS)
    def test_hgvs_Variant_equivalence_via_protein(self, hgvs: str, variant: Variant):
        """Test hgvs and Variant equivalence by comparing the protein prediction

        The goal here is to verify that the model representation of the Variant
        is usable by mutalyzer in the same way as the original HGVS description
        """
        d = _init_d(f"{self.transcript}:c.{hgvs}")
        variant_protein = protein_from_variant(variant, self.empty_transcript)
        description_protein = protein_from_description(d)

        assert variant_protein == description_protein

    VARIANTS = [
        ("10C>T", Variant(47576, 47577, inserted="A", deleted="G")),
        ("10del", Variant(47576, 47577, inserted="")),
        ("10_11insA", Variant(47576, 47576, inserted="T")),
        ("10_10delinsT", Variant(47576, 47577, inserted="A")),
        # Delete the original 'C', and insert CC (reverse complement)
        ("10dup", Variant(47576, 47577, inserted="GG")),
        ("10_10delinsCC", Variant(47576, 47577, inserted="GG")),
        # Delete the original "TC", and insert CTCT (reverse complement)
        ("10_11dup", Variant(47575, 47577, inserted="AGAG")),
        ("12_13delinsCTCT", Variant(47573, 47575, inserted="AGAG")),
        # Inversion, equivalent to 10C>G
        ("10_10inv", Variant(47576, 47577, inserted="C")),
        # Inversion, equivalent to 10_11delinsAG
        ("10_11inv", Variant(47575, 47577, inserted="CT")),
        # Inversion, not symetrical
        ("18_20inv", Variant(47566, 47569, inserted="GCA")),
        ("18_20delinsTGC", Variant(47566, 47569, inserted="GCA")),
        # Small mononucleotide repeat
        # 7_8T[4],  Variant representation = 7_8delinsTTTT
        ("7_8T[4]", Variant(47578, 47580, inserted="AAAA")),
        # Repeat of multiple distinct nucleotides
        # 10_13CT[4], Variant representation = 10_13delinsCTCTCTCT
        ("10_13CT[4]", Variant(47573, 47577, inserted="AGAGAGAG")),
    ]

    @pytest.mark.parametrize("hgvs, variant", VARIANTS)
    def test_hgvs_Variant_delins_model(self, hgvs: str, variant: Variant):
        """Test hgvs and Variant delins model equivalence directly

        The goal here is to verify that the Variant.from_model structure is
        always a simple delins model, even for duplications, repeats and
        inversions.
        """
        d = _init_d(f"{self.transcript}:c.{hgvs}")
        delins_model = delins_from_description(d)

        seq = sequence_from_description(d)
        variant_model = Variant.from_model(delins_model, sequence=seq)
        assert variant_model == variant


def manual(variant):
    tf = TestForward()
    tf = TestReverse()
    hgvs = f"{tf.transcript}:c.{variant}"
    d = _init_d(hgvs)

    delins_model = delins_from_description(d)
    print("=" * 10, "DELINS MODEL", "=" * 10)
    pprint(delins_model)

    seq = sequence_from_description(d)
    v = Variant.from_model(delins_model, sequence=seq)

    print("=" * 10, "VARIANT MODEL", "=" * 10)
    pprint(v.to_model())

    variant_protein = protein_from_variant(v, tf.empty_transcript)
    description_protein = protein_from_description(d)

    print(variant_protein, v)
    print(description_protein, d)
    print(v)
    print(variant_protein == description_protein)
    exit()
    print(protein_from_description(d))
    print(protein_from_variant(v, self.empty_transcript))
    exit()
    pprint(v.to_model())


if __name__ == "__main__":
    manual(sys.argv[1])
