"""
Module that contains the command line app, so we can still import __main__
without executing side effects
"""

from .ucsc import lookup_transcript
import argparse

import json


def logger_setup() -> None:
    import logging

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    logger.addHandler(ch)


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

    logger_setup()

    if args.transcript_id:
        r = lookup_transcript(args.transcript_id)
        print(r.json())


if __name__ == "__main__":
    main()
