#!/usr/bin/env python3

import itertools

from gtgt.mutalyzer import init_description, sequence_from_description

import sys

from gtgt.variant import Variant

sdhd = [
    "NM_003002.4",
    "NC_000011.10(NM_003002.4)",
    "ENST00000375549.8"
]

wt1 = [
    "NM_024426.6",
    "NC_000011.10(NM_024426.6)",
    "ENST00000452863.10"
]
data = {
    "sdhd": sdhd,
    "wt1": wt1
}

Range = tuple[int, int]

def _variant_size(variant: Variant) -> int:
    """
    Determine the net size of a variant
    """
    return int(variant.start - variant.end + len(variant.inserted))

def mutate_features(features: list[Range], variants: list[Variant]) -> list[Range]:
    """
    Mutate features based on the variants
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

            # If the variant is before the position, we increase the offset
            # for this position by the size of the variant
            if position > variant.end:
                offsets[index] += size

    # Apply the offsets
    boundaries = [position + offset for position, offset in zip(boundaries, offsets)]

    # Group into Ranges
    new_features = list()
    for i in range(0, len(boundaries), 2):
        new_features.append((boundaries[i], boundaries[i+1]))

    return new_features

def mutate_sequence(sequence: str, variants:list[Variant]) -> str:
    """Mutate the sequence with the variants"""
    if not variants:
        return sequence

    # Check that the variants are sorted
    prev = variants[0].end
    for variant in variants[1:]:
        if variant.start < prev:
            raise ValueError("Variants must be sorted")
        prev = variant.end

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
    "T": "A"
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
        hgvs = f"{transcript}:c.[10del; 16_17delinsT]"
        hgvs = f"{transcript}:c.10_11insA"
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
        print(variants)

        # print(d.get_selector_model()["exon"])
        # print(d.get_selector_model()["cds"])

        seq = get_sequence(sequence, d.get_selector_model()["exon"], d.is_inverted())
        print(f">{hgvs} {d.is_inverted()}")
        print(seq)

        print()
        
