import gzip
import re
from collections.abc import Iterator
from itertools import zip_longest
from typing import Any


def _parse_line(line: str) -> dict[str, Any]:
    """Parse a line from a uniprot selected.tab file"""
    header = "UniProtKB-AC UniProtKB-ID GeneID(EntrezGene) RefSeq GI PDB GO UniRef100 UniRef90 UniRef50 UniParc PIR NCBI-taxon MIM UniGene PubMed EMBL EMBL-CDS Ensembl Ensembl_TRS Ensembl_PRO Additional_PubMed".split(
        " "
    )
    # Fields that can hold multiple values
    multiple_values = "GeneID(EntrezGene) RefSeq GI PDB GO PIR MIM PubMed EMBL EMBL-CDS Ensembl Ensembl_TRS Ensembl_PRO Additional_PubMed".split()
    spline = line.strip("\n").split("\t")

    results = dict()
    for field, value in zip_longest(header, spline):
        if field in multiple_values:
            values = value.split("; ")
            # One or more values can be missing, this is indicated with a "-"
            values = [x if x != "-" else None for x in values]
        else:
            values = value
        # If the value is empty, store None instead
        results[field] = values if values else None
    return results


def _regex_reader(fname: str, regex: str) -> Iterator[str]:
    """Read lines that match the specified regex from a gzip file"""
    prog = re.compile(regex)
    with gzip.open(fname, "rt") as fin:
        for line in fin:
            if prog.search(line):
                yield line


def uniprot_id(fname: str, identifier: str) -> str | None:
    """Extract the uniprot ID for identifier from fname


    fname is expected to be HUMAN_9606_idmapping_selected.tab.gz, from
    https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

    identifier can be one of:
    - A RefSeq identifier (NP_, XP_, etc)
    - An Ensembl transcript id (ENST)
    """
    for line in _regex_reader(fname, identifier):
        data = _parse_line(line)
        if identifier in data["Ensembl_TRS"] or identifier in data["RefSeq"]:
            return str(data["UniProtKB-AC"])
    else:
        return None
