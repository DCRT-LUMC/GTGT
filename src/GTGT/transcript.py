from .bed import Bed

from typing import List


class Transcript:
    def __init__(self, records: List[Bed]):
        for record in records:
            if record.name == "exons":
                self.exons = record
            elif record.name == "cds":
                self.cds = record
