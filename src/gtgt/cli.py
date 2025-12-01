"""
Module that contains the command line app, so we can still import __main__
without executing side effects
"""

import argparse
import dataclasses
import json
import logging
import os
import secrets

from gtgt.flask import render
from gtgt.models import TranscriptModel
from gtgt.transcript import Result, Transcript

from .mutalyzer import (
    Variant,
    generate_therapies,
    init_description,
    sequence_from_description,
)
from .provider import Provider
from .variant_validator import lookup_variant


def set_logging(level: str) -> None:
    format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    format = "%(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(
        format=format,
        level=level.upper(),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Description of command.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--log",
        choices=["debug", "info", "warning", "error", "critical"],
        default="info",
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

    export_parser = subparsers.add_parser(
        "export", help="Export the specified Transcript to BED format"
    )
    export_parser.add_argument(
        "hgvs", help="HGVS description of the transcript of interest"
    )

    render_parser = subparsers.add_parser(
        "render", help="Render the HTML template (helper)"
    )

    render_parser.add_argument("--results", type=str, default="")
    render_parser.add_argument("--error", "-e", help="Error payload")
    render_parser.add_argument("--variant", type=str, default="10A>T")
    args = parser.parse_args()

    set_logging(args.log)
    logger = logging.getLogger(__name__)

    provider = Provider(args.cachedir)

    if args.command == "transcript":
        d = init_description(f"{args.transcript_id}:c.=")
        t = Transcript.from_description(d)
        ts = TranscriptModel.from_transcript(t)
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
            logger.critical(
                "Missing modules, please install with 'pip install gtgt[api_server]'"
            )
            exit(-1)
        uvicorn.run(app, host=args.host)
    elif args.command == "mutate":
        desc = f"{args.transcript_id}:c.="
        d = init_description(desc)
        for therapy in generate_therapies(d):
            print(f"{therapy.name}: {therapy.hgvsc}")
    elif args.command == "analyze":
        d = init_description(args.hgvs)
        transcript = Transcript.from_description(d)
        # Convert Result objects to dict
        results = [dataclasses.asdict(x) for x in transcript.analyze(args.hgvs)]
        print(json.dumps(results, indent=True, default=vars))
    elif args.command == "webserver":
        try:
            from .flask import app as flask_app
        except ModuleNotFoundError:
            logger.critical(
                f"Missing modules, please install with 'pip install gtgt[webserver]'"
            )
            exit(-1)
        if not flask_app.config.get("SECRET_KEY"):
            flask_app.secret_key = secrets.token_hex()
        flask_app.run(args.host, debug=args.debug)
    elif args.command == "export":
        # Get the transcript
        transcript = Transcript.from_description(init_description(args.hgvs))

        # Mutate the transcript
        d = init_description(args.hgvs)
        sequence = sequence_from_description(d)
        input_variants = [
            Variant.from_model(delins, sequence=sequence)
            for delins in d.delins_model["variants"]
        ]
        transcript.mutate(d, input_variants)

        for record in transcript.records():
            print(record)
    elif args.command == "render":
        if args.results:
            with open(args.results) as fin:
                file_payload = [Result.from_dict(x) for x in json.load(fin)]
        else:
            file_payload = []

        template_file = "templates/index.html.j2"
        print(
            render(
                template_file,
                variant=args.variant,
                results=file_payload,
                error=args.error,
            )
        )
    else:
        parser.print_help()
        exit(1)


if __name__ == "__main__":
    main()
