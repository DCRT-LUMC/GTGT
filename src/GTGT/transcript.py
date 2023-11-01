from copy import deepcopy

from .bed import Bed

from typing import List


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
        for record in self.records():
            record.intersect(selector)

    def overlap(self, selector: Bed) -> None:
        for record in self.records():
            record.overlap(selector)

    def subtract(self, selector: Bed) -> None:
        pass
