#!/usr/bin/env python3

import pytest
import sys
from gtgt.mutalyzer import Variant
from mutalyzer.description import Description
from mutalyzer.protein import get_protein_description
from mutalyzer.reference import get_protein_selector_model
from mutalyzer.description_model import get_reference_id, variants_to_description
import json


def pprint(thing):
    print(json.dumps(thing, indent=True))

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
    """ Return the internal delins model from a Description object """
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

    desc, ref, predicted, *rest = get_protein_description(delins_models, d.references, selector_model)
    return predicted

class TestForward():
    """Test equivalent models on the forward strand"""
    # SDHD, forward strand
    transcript="ENST00000375549.8"
    empty_transcript=_init_d(f"{transcript}:c.=")

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
    ]

    @pytest.mark.parametrize("hgvs,variant", VARIANTS)
    def test_hgvs_Variant_equivalence_via_protein(self, hgvs: str, variant: Variant):
        """Test hgvs and Variant equivalence by comparing the protein prediction"""
        d = _init_d(f"{self.transcript}:c.{hgvs}")
        variant_protein = protein_from_variant(variant, self.empty_transcript)
        description_protein = protein_from_description(d)

        assert variant_protein == description_protein

def manual(variant):
    tf = TestReverse()
    hgvs = f"{tf.transcript}:c.{variant}"
    d = _init_d(hgvs)

    delins_model = delins_from_description(d)

    seq = sequence_from_description(d)
    v = Variant.from_model(delins_model, sequence=seq)

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
