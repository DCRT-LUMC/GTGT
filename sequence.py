#!/usr/bin/env python3

from copy import copy
import itertools
from typing import Iterable, Sequence

from gtgt.mutalyzer import init_description, sequence_from_description

import sys

from gtgt.range import intersect
from gtgt.variant import Variant
from Bio.Seq import Seq

import time
import resource

sdhd = [
    "NC_000011.10(NM_003002.4)",
    "NM_003002.4",
    "ENST00000375549.8"
]

wt1 = [
    "NC_000011.10(NM_024426.6)",
    "NM_024426.6",
    "ENST00000452863.10"
]
data = {
    "sdhd": sdhd,
    "wt1": wt1
}

Range = tuple[int, int]

class MiniTranscript:
    """
    Minimalistic transcript

    If you make a copy of MiniTranscript objects, make sure to replace features
    and not modify them, since members are shared between copies by default.
    """
    def __init__(self, inverted: bool, sequence:str, exons: list[Range], cds: list[Range] | None = None, variants: list[Variant] | None = None) -> None:
        self.inverted = inverted
        self.sequence = sequence
        self.exons = exons
        self.cds = cds
        self.variants = tuple(variants) if variants else tuple()

    def get_sequence(self) -> str:
        """Return the sequence of the transcript"""
        return self.sequence

    def _get_features_sequence(self, sequence:str, features: list[Range]) -> str:
        """Return the sequence of the features, in forward orientation"""
        segments = list()
        for start, end in features:
            segments.append(sequence[start:end])

        if self.inverted:
            return reverse_complement(''.join(segments))
        else:
            return ''.join(segments)

    def get_exon_sequence(self) -> str:
        """Return the exon sequence in forward orientation"""
        exons = mutate_features(self.exons, self.variants)
        sequence = mutate_sequence(self.sequence, self.variants)

        return self._get_features_sequence(sequence, exons)

    def get_coding_sequence(self) -> str:
        """Return the coding sequence in forward orientation"""
        if not self.cds:
            raise ValueError("Transcript does not have a cds")
        exons = mutate_features(self.exons, self.variants)
        cds = mutate_features(self.cds, self.variants)[0]
        sequence = mutate_sequence(self.sequence, self.variants)

        coding = list()
        for exon in exons:
            coding += intersect(exon, cds)

        return self._get_features_sequence(sequence, coding)

    def get_protein_sequence(self) -> str:
        """Return the protein sequence for the transcript"""
        coding_sequence = self.get_coding_sequence()
        # Add padding if the coding sequence is not divisible by three
        if (padding:=3 - len(coding_sequence)%3):
            coding_sequence += "N"*padding

        predicted = str(Seq(coding_sequence).translate()) # type: ignore
        first_stop = predicted.find("*")
        return predicted[: first_stop + 1]



def _variant_size(variant: Variant) -> int:
    """
    Determine the net size of a variant
    """
    return int(variant.start - variant.end + len(variant.inserted))

def mutate_features(features: Sequence[Range], variants: Sequence[Variant]) -> list[Range]:
    """
    Mutate features based on the variants

    The variants do not have to be sorted
    """
    # Flatten the features, since we can act on them individually
    boundaries = list(itertools.chain.from_iterable(features))

    # Keep a running tally of the offsets for each boundary, because we cannot
    # update the boundaries until we have seen all the variants
    offsets = [0 for _ in boundaries]

    for variant in variants:
        size = _variant_size(variant)
        for index in range(len(boundaries)):
            position = boundaries[index]

            # If the variant is before or at the position, we increase the offset
            # for this position by the size of the variant
            if position >= variant.end:
                offsets[index] += size

    # Apply the offsets
    boundaries = [position + offset for position, offset in zip(boundaries, offsets)]

    # Group into Ranges
    new_features = list()
    for i in range(0, len(boundaries), 2):
        new_features.append((boundaries[i], boundaries[i+1]))

    return new_features

def mutate_sequence(sequence: str, variants:Sequence[Variant]) -> str:
    """
    Mutate sequence based on the variants

    Note: the variants must be sorted
    """
    if not variants:
        return sequence

    # Check that the variants are sorted
    prev = variants[0].end
    for variant in variants[1:]:
        if variant.start < prev:
            raise ValueError("Variants must be sorted")
        prev = variant.end

    # The end of the last variant cannot be beyond the end of the sequence
    if variants[-1].end > len(sequence):
        raise ValueError(f"{variants[-1]} lies beyond the end of the sequence")

    # Store the segments of the mutated sequence
    segments = list()

    # Keep track of the last position in the sequence we have included
    position = 0

    for variant in variants:
        # Store the segment before the current variant
        segments.append(sequence[position:variant.start])
        segments.append(variant.inserted)
        position = variant.end
    else:
        segments.append(sequence[position:])

    return "".join(segments)

mapping = str.maketrans({
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "a": "u",
    "c": "g",
    "g": "c",
    "u": "a",
})
def reverse_complement(seq: str) -> str:
    return seq.translate(mapping)[::-1]

def get_sequence(seq: str, blocks: list[Range], reverse: bool) -> str:
    """
    Get sequence
    """
    new = ''.join((seq[start:end] for start, end in blocks))

    if reverse:
        return reverse_complement(new)
    else:
        return new

if __name__ == "__main__":
    gene = sys.argv[1]

    for transcript in data[gene]:
        # Make DNA
        hgvs = f"{transcript}:r.[10del; 16_17delinsT]"
        # hgvs = f"{transcript}:c.10_11insA"
        # hgvs = f"{transcript}:r.[10del; 20A>T]"
        d = init_description(hgvs)
        # print(f"{hgvs}")
        # print(f"{d.is_inverted()=}")

        # Get the Variants on the internal coordinate system
        sequence = sequence_from_description(d)
        variants = [
            Variant.from_model(delins, sequence=sequence)
            for delins in d.delins_model["variants"]
        ]

        exons=d.get_selector_model()["exon"]
        cds=d.get_selector_model()['cds']

        print(variants)
        print(exons, cds)

        t = MiniTranscript(d.is_inverted(), sequence, exons, cds=cds, variants=variants)

        # Make a wildtype version of the transcript
        wt = copy(t)
        wt.variants = ()

        # print(d.get_selector_model()["exon"])
        # print(d.get_selector_model()["cds"])

        seq = get_sequence(sequence, d.get_selector_model()["exon"], d.is_inverted())
        print(f">{hgvs} {d.is_inverted()}")
        print('mutalyzer ', d.protein['predicted'])
        print('mut WT    ', d.protein['reference'])
        print('minitrans ', t.get_protein_sequence())
        print('mini WT   ', wt.get_protein_sequence())
        print(t.sequence)
        print()

        print("Before")
        print(f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss:,} KiB")
        # time.sleep(10)

        # Make 10 copies of the transcript
        copies = [copy(t) for _ in range(10)]
        print("After")
        print(f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss:,} KiB")
        
