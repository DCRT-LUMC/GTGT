#!/usr/bin/env python3

from copy import copy
import itertools
from typing import Iterable, Sequence

from mutalyzer_crossmapper import crossmapper

from gtgt.mutalyzer import init_description, sequence_from_description
from gtgt.transcript import MiniTranscript, reverse_complement

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
        hgvs="ENST00000436367.6:c.-730+100A>T"
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
        mapper = crossmapper.Coding(exons, cds[0], inverted=d.is_inverted())

        for var in variants:
            print(var)
            s = mapper.coordinate_to_coding(var.start)
            e = mapper.coordinate_to_coding(var.end)
            print(f"{s}_{e}{var.inserted}")
        exit()

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
        
