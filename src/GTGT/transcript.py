from copy import deepcopy

from .bed import Bed

from typing import List, Dict, Tuple

class Transcript:
    def __init__(self, exons: Bed, cds: Bed):
        self.exons = exons
        self.cds = cds

        # Set the coding region
        coding = deepcopy(self.exons)
        coding.name = "coding"
        coding.intersect(self.cds)
        self.coding = coding

    def records(self) -> List[Bed]:
        """Return the Bed records that make up the Transcript"""
        return [self.exons, self.cds, self.coding]

    def intersect(self, selector: Bed) -> None:
        """Update transcript to only contain features that intersect the selector"""
        for record in self.records():
            record.intersect(selector)

    def overlap(self, selector: Bed) -> None:
        """Update transcript to only contain features that overlap the selector"""
        for record in self.records():
            record.overlap(selector)

    def subtract(self, selector: Bed) -> None:
        """Remove all features from transcript that intersect the selector"""
        for record in self.records():
            record.subtract(selector)

    def exon_skip(self, selector: Bed) -> None:
        """Remove the exon(s) that overlap the selector from the transcript"""
        exons_to_skip = deepcopy(self.exons)
        exons_to_skip.overlap(selector)
        self.subtract(exons_to_skip)

    def compare(self, other: object) -> Dict[str, float]:
        """Compare the size of each record in the transcripts"""
        if not isinstance(other, Transcript):
            raise NotImplementedError

        # Compare each record that makes up self and other
        # The comparison will fail if the record.name does not match
        cmp = dict()
        for record1, record2 in zip(self.records(), other.records()):
            cmp[record1.name] = record1.compare(record2)

        return cmp


    def _reverse_cdot_to_genomic(self, cdot: str) -> int:
        """
        Convert a HGVS c. location to genomic, for transcripts on the
        reverse strand
        """
        pos = int(cdot) - 1
        return self.cds.chromEnd - pos

    def cdot_to_genomic(self, cdot: str) -> int:
        """Convert a HGVS c. location to genomic"""
        if self.exons.strand == "-":
            return self._reverse_cdot_to_genomic(cdot)
        elif self.exons.strand == "+":
            pass
        else:
            raise RuntimeError("Unable to convert c. position, transcript has no strand")
        return 0

def find_distance_endpoint(ranges: List[range], start: int, distance: int) -> int:
    """
    Travel over a list of ranges and determine the endpoint

    Starting from the 'start' position, travel 'distance' over a list of ranges
    and report the endpoint
    """
    # Find the start position in the ranges
    for index, range_ in enumerate(ranges):
        if start in range_:
            break
    else:
        raise ValueError(f"Start position '{start}' lies outside the ranges")

    # check if we can travel distance in the current range
    endpoint = start + distance
    if endpoint in range_:
        return ranges[index].start + start + distance
    else:
        new_start = range_.start
        new_distance = None
    # If not, check which way we go over te boundary
    # Update index to match the new range
    # Recurse
    return 0
