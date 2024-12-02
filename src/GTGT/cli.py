"""
Module that contains the command line app, so we can still import __main__
without executing side effects
"""

from .wrappers import lookup_transcript
from .variant_validator import lookup_variant
from .provider import Provider
from .mutalyzer import exonskip
from .mutalyzer import HGVS

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
    subparsers = parser.add_subparsers(dest="command")
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

    mutator_parser = subparsers.add_parser(
        "mutate", help="Mutate the specified transcript"
    )

    mutator_parser.add_argument("transcript_id", help="The transcript to mutate")

    args = parser.parse_args()

    logger = logger_setup()

    provider = Provider(args.cachedir)

    if args.command == "transcript":
        ts = lookup_transcript(provider, args.transcript_id)
        print(ts.model_dump_json())
    elif args.command == "links":
        logger.debug(args)
        links = lookup_variant(provider, args.hgvs_variant).url_dict()
        for website, url in links.items():
            print(f"{website}: {url}")
    elif args.command == "server":
        try:
            from .app import app, uvicorn
        except ModuleNotFoundError:
            print("Missing modules, please install with 'pip install gtgt[server]'")
            exit(-1)
        uvicorn.run(app, host=args.host)
    elif args.command == "mutate":
        desc = f"{args.transcript_id}:c.="
        hgvs_transcript = HGVS(description=desc)
        for hgvs in exonskip(hgvs_transcript):
            print(hgvs.description)
        print(args.transcript_id)
    else:
        raise NotImplementedError


if __name__ == "__main__":
    main()
