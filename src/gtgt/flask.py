import os
from collections.abc import Sequence
from typing import Any, Mapping

import mutalyzer_hgvs_parser
from flask import Flask
from jinja2 import Environment, FileSystemLoader

from gtgt.mutalyzer import init_description
from gtgt.transcript import Result

from . import Transcript
from .provider import Provider
from .variant_validator import lookup_variant

hgvs_error = (
    mutalyzer_hgvs_parser.exceptions.UnexpectedCharacter,
    mutalyzer_hgvs_parser.exceptions.UnexpectedEnd,
)

app = Flask(__name__)
provider = Provider()


payload = dict[str, Any] | None


def render(
    template_file: str,
    variant: str | None = None,
    results: Sequence[Result] = [],
    links: Mapping[str, str] | None = None,
    error: Mapping[str, str] | None = None,
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
        d["results"] = organize_results(results)
    if links:
        d["links"] = links
    if error:
        d["error"] = error

    return template.render(**d)


def organize_results(
    results: Sequence[Result],
) -> dict[str, list[Result] | Result]:
    """
    The Results are grouped in a dict with the following keys
    - input, which contains the Result for the input
    - modified, which contains the Results where one or more of the input
      varianst have been changed
    - all, which contains all other Results, including the wildtype etc, in
      the original order
    """

    # Keep track of the Results we still have to precess
    todo = [x for x in results]

    # Extract the input and wildtype results, wich are special
    input_ = next((x for x in results if x.therapy.name == "Input"))
    wildtype = next((x for x in results if x.therapy.name == "Wildtype"))
    todo.remove(input_)
    todo.remove(wildtype)

    # Extract the variants from the input
    input_variants = input_.therapy.variants

    # Find all Results where one or more of the input variants are modified or removed
    modified = list()
    for result in todo:
        for variant in input_variants:
            if variant not in result.therapy.variants:
                modified.append(result)
                break

    # Exclude the Results that modify the input variant
    todo = [x for x in todo if x not in modified]

    return {
        "input": [input_, wildtype],
        "modified": modified,
        "all": todo,
    }


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

    return error


@app.route("/")
@app.route("/<variant>")
def result(variant: str | None = None) -> str:
    template_file = "templates/index.html.j2"

    # If no variant was specified
    if not variant:
        return render(template_file)

    # Invalid user input
    if error := validate_user_input(variant):
        return render(template_file, variant=variant, error=error)

    # Analyze the transcript
    try:
        transcript = Transcript.from_description(init_description(variant))
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
