from .variant_validator import lookup_variant
from .provider import Provider
from .wrappers import lookup_transcript

from flask import Flask, render_template
from typing import Optional, Dict

import mutalyzer_hgvs_parser

hgvs_error = (
    mutalyzer_hgvs_parser.exceptions.UnexpectedCharacter,
    mutalyzer_hgvs_parser.exceptions.UnexpectedEnd,
)

app = Flask(__name__)
provider = Provider()


def validate_user_input(input: str) -> Dict[str, str]:
    """
    Validate the user input

    If there is an error, return a dict with summary and details of the error
    """
    error = dict()

    # Test if the variant is valid HGVS
    try:
        model = mutalyzer_hgvs_parser.to_model(input)
    except hgvs_error as e:
        error["summary"] = "Not a valid HGVS description"
        error["details"] = str(e)
        return error

    if not input.startswith("ENST"):
        error["summary"] = "Not an ensembl transcript"
        error["details"] = "Currently, only ensembl transcripts (ENST) are supported"
        return error

    if model["coordinate_system"] != "c":
        error["summary"] = "Only 'c.' coordinates are supported"
        error["details"] = "Please convert your variant to the 'c.' coordinate system"
        return error

    return error


@app.route("/")
@app.route("/<variant>")
def result(variant: Optional[str] = None) -> str:
    template_file = "index.html.j2"

    # If no variant was specified
    if not variant:
        return render_template(template_file)

    # Invalid user input
    if error := validate_user_input(variant):
        return render_template(template_file, variant=variant, error=error)

    # Analyze the transcript
    try:
        transcript_id = variant.split(":")[0]
        transcript_model = lookup_transcript(provider, transcript_id)
        transcript = transcript_model.to_transcript()
        results = transcript.analyze(variant)
    except Exception as e:
        error = {"summary": "Analysis failed", "details": str(e)}
        results = None

    # Get external links
    try:
        links = lookup_variant(provider, variant).url_dict()
    except Exception as e:
        links = dict()

    if error:
        return render_template(template_file, variant=variant, error=error)
    else:
        return render_template(
            template_file, results=results, links=links, variant=variant
        )
