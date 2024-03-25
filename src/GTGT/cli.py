"""
Module that contains the command line app, so we can still import __main__
without executing side effects
"""

from .ensembl import lookup_transcript
from .ucsc import lookup_knownGene
from .bed import Bed
from .provider import Provider

import argparse
import json
import logging
from .models import BedModel


def logger_setup() -> logging.Logger:
    import logging

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    logger.addHandler(ch)

    return logger


def main() -> None:
    parser = argparse.ArgumentParser(description="Description of command.")
    subparsers = parser.add_subparsers()
    parser.add_argument("--name", default="world", required=False)

    transcript_parser = subparsers.add_parser(
        "transcript", help="Transcript Information"
    )

    transcript_parser.add_argument(
        "transcript_id", type=str, help="Transcript of interest"
    )
    args = parser.parse_args()

    logger = logger_setup()

    if args.transcript_id:
        provider = Provider()
        r = lookup_transcript(provider, args.transcript_id)
        logger.debug(r)
        track = lookup_knownGene(provider, r)
        print(r.json())
        knownGene = track["knownGene"][0]
        print(json.dumps(knownGene))
        bm = BedModel.from_ucsc(knownGene)
        print(bm.model_dump_json())
        exit()
        # bed = Bed.from_ucsc(knownGene)
        # print(bed)
        exit()

        print(json.dumps(track))


if __name__ == "__main__":
    main()
