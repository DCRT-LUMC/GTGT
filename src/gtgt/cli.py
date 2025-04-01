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

    api_server_parser = subparsers.add_parser(
        "api_server", help="Run the GTGT API server"
    )
    api_server_parser.add_argument(
        "--host", default="localhost", help="Hostname to listen on"
    )

    web_server_parser = subparsers.add_parser(
        "webserver", help="Run the GTGT web server"
    )
    web_server_parser.add_argument(
        "--host", default="localhost", help="Hostname to listen on"
    )
    web_server_parser.add_argument(
        "--debug", default=False, action="store_true", help="Run Flask in debug mode"
    )

    mutator_parser = subparsers.add_parser(
        "mutate", help="Mutate the specified transcript"
    )

    mutator_parser.add_argument("transcript_id", help="The transcript to mutate")

    analyze_parser = subparsers.add_parser(
        "analyze", help="Analyze all possible exons skips for the spcified HGVS variant"
    )
    analyze_parser.add_argument(
        "hgvs", help="HGVS description of the transcript of interest"
    )
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
    elif args.command == "api_server":
        try:
            from .app import app, uvicorn
        except ModuleNotFoundError:
            print("Missing modules, please install with 'pip install gtgt[api_server]'")
            exit(-1)
        uvicorn.run(app, host=args.host)
    elif args.command == "mutate":
        desc = f"{args.transcript_id}:c.="
        hgvs_transcript = HGVS(description=desc)
        for hgvs in exonskip(hgvs_transcript):
            print(hgvs.description)
        print(args.transcript_id)
    elif args.command == "analyze":
        transcript_id = args.hgvs.split(":")[0]
        transcript_model = lookup_transcript(provider, transcript_id)
        transcript = transcript_model.to_transcript()
        results = transcript.analyze(args.hgvs)
        print(json.dumps(results, indent=True))
    elif args.command == "webserver":
        try:
            from .flask import app as flask_app
        except ModuleNotFoundError:
            print("Missing modules, please install with 'pip install gtgt[webserver]'")
            exit(-1)
        flask_app.run(args.host, debug=args.debug)
    else:
        raise NotImplementedError()


if __name__ == "__main__":
    main()
