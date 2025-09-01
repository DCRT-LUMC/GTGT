from collections.abc import Sequence
import os
from typing import Any, Mapping, Optional
from gtgt.transcript import Result

import mutalyzer_hgvs_parser
from flask import Flask
from jinja2 import Environment, FileSystemLoader

from .provider import Provider
from .variant_validator import lookup_variant
from .wrappers import lookup_transcript

hgvs_error = (
    mutalyzer_hgvs_parser.exceptions.UnexpectedCharacter,
    mutalyzer_hgvs_parser.exceptions.UnexpectedEnd,
)

app = Flask(__name__)
provider = Provider()


payload = Optional[dict[str, Any]]


def render(
    template_file: str,
    variant: Optional[str] = None,
    results: Sequence[Result] = [],
    links: Optional[Mapping[str, str]] = None,
    error: Optional[Mapping[str, str]] = None,
) -> str:
    """Render the html template using jinja2 directly"""
    root = os.path.dirname(__file__)
    template_file = os.path.join(root, template_file)
    template_folder = os.path.dirname(template_file)
    environment = Environment(loader=FileSystemLoader(template_folder))
    template = environment.get_template(os.path.basename(template_file))

    # Exclude any payloads that are not set
    d: dict[str, Any] = dict()
    if variant:
        d["variant"] = variant
    if results:
        d["results"] = results
    if links:
        d["links"] = links
    if error:
        d["error"] = error

    return template.render(**d)


def validate_user_input(input: str) -> Mapping[str, str]:
    """
    Validate the user input

    If there is an error, return a dict with summary and details of the error
    """
    error = dict()

    # Test if the variant is valid HGVS
    try:
        mutalyzer_hgvs_parser.to_model(input)
    except hgvs_error as e:
        error["summary"] = "Not a valid HGVS description"
        error["details"] = str(e)
        return error

    if not input.startswith("ENST"):
        error["summary"] = "Not an ensembl transcript"
        error["details"] = "Currently, only ensembl transcripts (ENST) are supported"
        return error

    return error


@app.route("/")
@app.route("/<variant>")
def result(variant: Optional[str] = None) -> str:
    template_file = "templates/index.html.j2"

    # If no variant was specified
    if not variant:
        return render(template_file)

    # Invalid user input
    if error := validate_user_input(variant):
        return render(template_file, variant=variant, error=error)

    # Analyze the transcript
    try:
        transcript_id = variant.split(":")[0]
        transcript_model = lookup_transcript(provider, transcript_id)
        transcript = transcript_model.to_transcript()
        results = transcript.analyze(variant)
    except Exception as e:
        error = {"summary": "Analysis failed", "details": str(e)}
        results = []

    # Get external links
    try:
        links = lookup_variant(provider, variant).url_dict()
    except Exception as e:
        links = dict()

    if error:
        return render(template_file, variant=variant, error=error)
    else:
        return render(template_file, results=results, links=links, variant=variant)
