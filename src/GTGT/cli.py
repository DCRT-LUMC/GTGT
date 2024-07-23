"""
Module that contains the command line app, so we can still import __main__
without executing side effects
"""

from .ensembl import lookup_transcript
from .ucsc import lookup_knownGene
from .bed import Bed
from .variant_validator import lookup_variant
from .provider import Provider
from .models import BedModel, TranscriptModel

import argparse
import json
import logging
import os


def logger_setup() -> logging.Logger:
    import logging

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    logger.addHandler(ch)

    return logger


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Description of command.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers()
    parser.add_argument("--cachedir", type=str, default=os.environ.get("GTGT_CACHE"))

    transcript_parser = subparsers.add_parser(
        "transcript", help="Transcript Information"
    )

    transcript_parser.add_argument(
        "transcript_id", type=str, help="Transcript of interest"
    )

    link_parser = subparsers.add_parser("links", help="Links to external resources")
    link_parser.add_argument("hgvs_variant", type=str, help="Variant of interest")

    server_parser = subparsers.add_parser("server", help="Run the GTGT server")
    server_parser.add_argument(
        "--host", default="0.0.0.0", help="Hostname to listen on"
    )

    args = parser.parse_args()

    logger = logger_setup()

    provider = Provider(args.cachedir)

    if "transcript_id" in args:
        r = lookup_transcript(provider, args.transcript_id)
        logger.debug(r)
        track_name = "ncbiRefSeq"
        track_name = "ensGene"
        track = lookup_knownGene(provider, r, track_name)
        knownGene = track[track_name][0]
        bm = BedModel.from_ucsc(knownGene)

        exons = bm.copy()
        cds = BedModel.from_bed(Bed(bm.chrom, bm.thickStart, bm.thickEnd))
        ts = TranscriptModel(exons=exons, cds=cds)
        print(ts.model_dump_json())
    elif "hgvs_variant" in args:
        logger.debug(args)
        links = lookup_variant(provider, args.hgvs_variant)
        for url in links.url("omim"):
            print("OMIM:", url)
        print("LOVD:", links.url("lovd"))
        print("GTEx:", links.url("gtex"))
        print("UniProt", links.url("uniprot"))
        print("DECIPHER", links.url("decipher"))
        print("ClinVar", links.url("clinvar"))
        print("HGNC", links.url("hgnc"))
        print("UCSC", links.url("ucsc"))
    elif "host" in args:
        try:
            from .app import app, uvicorn
        except ModuleNotFoundError:
            print("Missing modules, please install with 'pip install gtgt[server]'")
            exit(-1)
        uvicorn.run(app, host=args.host)
    else:
        raise NotImplementedError


if __name__ == "__main__":
    main()
