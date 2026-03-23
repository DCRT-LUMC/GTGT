"""
Module that contains the command line app, so we can still import __main__
without executing side effects
"""

import argparse
import dataclasses
import json
import logging
import secrets
from typing import Any

from gtgt.flask import render as flask_render
from gtgt.transcript import Result, Transcript

from .mutalyzer import (
    init_description,
    sequence_from_description,
)
from .therapy import generate_therapies
from .variant import Variant
from .variant_validator import lookup_variant

logger = logging.getLogger(__name__)


def set_logging(level: str) -> None:
    format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    format = "%(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(
        format=format,
        level=level.upper(),
    )


def _set_transcript_parser(subparsers: Any) -> None:
    transcript_parser = subparsers.add_parser(
        "transcript", help="Transcript Information"
    )

    transcript_parser.add_argument(
        "transcript_id", type=str, help="Transcript of interest"
    )
    transcript_parser.set_defaults(func=transcript)


def _set_link_parser(subparsers: Any) -> None:
    link_parser = subparsers.add_parser("links", help="Links to external resources")
    link_parser.add_argument("hgvs_variant", type=str, help="Variant of interest")
    link_parser.set_defaults(func=link)


def _set_webserver_parser(subparsers: Any) -> None:
    web_server_parser = subparsers.add_parser(
        "webserver", help="Run the GTGT web server"
    )
    web_server_parser.add_argument(
        "--host", default="localhost", help="Hostname to listen on"
    )
    web_server_parser.add_argument(
        "--debug", default=False, action="store_true", help="Run Flask in debug mode"
    )
    web_server_parser.set_defaults(func=webserver)


def _set_mutator_parser(subparsers: Any) -> None:
    mutator_parser = subparsers.add_parser(
        "mutate", help="Mutate the specified transcript"
    )

    mutator_parser.add_argument("transcript_id", help="The transcript to mutate")
    mutator_parser.set_defaults(func=mutate)


def _set_analyze_parser(subparsers: Any) -> None:
    analyze_parser = subparsers.add_parser(
        "analyze", help="Analyze all possible exons skips for the spcified HGVS variant"
    )
    analyze_parser.add_argument(
        "hgvs", help="HGVS description of the transcript of interest"
    )
    analyze_parser.add_argument(
        "--protein-domains",
        help="Fetch supported protein domains from UCSC",
        default=False,
        action="store_true",
    )
    analyze_parser.add_argument(
        "--extended",
        help="Output all evaluated therapies",
        default=False,
        action="store_true",
    )
    analyze_parser.set_defaults(func=analyze)


def _set_export_parser(subparsers: Any) -> None:
    export_parser = subparsers.add_parser(
        "export", help="Export the specified Transcript to BED format"
    )
    export_parser.add_argument(
        "hgvs", help="HGVS description of the transcript of interest"
    )
    export_parser.add_argument(
        "--protein-domains",
        help="Fetch supported protein domains from UCSC",
        default=False,
        action="store_true",
    )
    export_parser.set_defaults(func=export)


def _set_render_parser(subparsers: Any) -> None:
    render_parser = subparsers.add_parser(
        "render", help="Render the HTML template (helper)"
    )

    render_parser.add_argument("--results", type=str, default="")
    render_parser.add_argument("--error", "-e", help="Error payload")
    render_parser.add_argument("--variant", type=str, default="10A>T")
    render_parser.set_defaults(func=render)


def transcript(args: argparse.Namespace) -> None:
    d = init_description(f"{args.transcript_id}:c.=")
    t = Transcript.from_description(d)
    print(t)


def link(args: argparse.Namespace) -> None:
    links = lookup_variant(args.hgvs_variant).url_dict()
    for website, url in links.items():
        print(f"{website}: {url}")


def mutate(args: argparse.Namespace) -> None:
    desc = f"{args.transcript_id}:c.="
    d = init_description(desc)
    for therapy in generate_therapies(d):
        print(f"{therapy.name}: {therapy.hgvsc}")


def analyze(args: argparse.Namespace) -> None:
    d = init_description(args.hgvs)
    transcript = Transcript.from_description(d)
    if args.protein_domains:
        transcript.lookup_protein_domains(d)
    # Convert Result objects to dict
    results = [
        dataclasses.asdict(x)
        for x in transcript.analyze(args.hgvs, extended=args.extended)
    ]
    print(json.dumps(results, indent=True, default=vars))


def webserver(args: argparse.Namespace) -> None:
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


def export(args: argparse.Namespace) -> None:
    # Get the transcript
    d = init_description(args.hgvs)
    transcript = Transcript.from_description(d)
    if args.protein_domains:
        transcript.lookup_protein_domains(d)

    # Mutate the transcript
    sequence = sequence_from_description(d)
    input_variants = [
        Variant.from_model(delins, sequence=sequence)
        for delins in d.delins_model["variants"]
    ]
    transcript.mutate(d, input_variants)

    for record in transcript.records():
        print(record)


def render(args: argparse.Namespace) -> None:
    if args.results:
        with open(args.results) as fin:
            file_payload = [Result.from_dict(x) for x in json.load(fin)]
    else:
        file_payload = []

    template_file = "templates/index.html.j2"
    print(
        flask_render(
            template_file,
            variant=args.variant,
            results=file_payload,
            error=args.error,
        )
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

    _set_transcript_parser(subparsers)
    _set_link_parser(subparsers)
    _set_webserver_parser(subparsers)
    _set_mutator_parser(subparsers)
    _set_analyze_parser(subparsers)
    _set_export_parser(subparsers)
    _set_render_parser(subparsers)

    args = parser.parse_args()

    set_logging(args.log)

    if args.command is None:
        parser.print_help()
        exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
