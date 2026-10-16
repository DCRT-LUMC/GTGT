import os
from collections.abc import Sequence
from typing import Any, Mapping

import mutalyzer_hgvs_parser
from flask import Flask, request
from jinja2 import Environment, FileSystemLoader

from gtgt.mutalyzer import init_description
from gtgt.transcript import Result

from . import Allele
from .variant_validator import lookup_variant

hgvs_error = (
    mutalyzer_hgvs_parser.exceptions.UnexpectedCharacter,
    mutalyzer_hgvs_parser.exceptions.UnexpectedEnd,
)

app = Flask(__name__)


payload = dict[str, Any] | None


def render(
    template_file: str,
    **kwargs: Any,
) -> str:
    """Render the html template using jinja2 directly"""
    root = os.path.dirname(__file__)
    template_file = os.path.join(root, template_file)
    template_folder = os.path.dirname(template_file)
    environment = Environment(loader=FileSystemLoader(template_folder))
    template = environment.get_template(os.path.basename(template_file))

    # Exclude any payloads that are not set
    if "results" in kwargs:
        kwargs["results"] = organize_results(kwargs["results"])

    return template.render(**kwargs)


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
def result() -> str:
    template_file = "templates/index.html.j2"

    # Extract parameter from the request
    protein_domains = request.args.get("protein_domains") == "true"
    extended = request.args.get("extended") == "true"

    dna_hgvs = request.args.get("dna")
    rna_hgvs = request.args.get("rna")

    # If no variant was specified
    if not dna_hgvs or not rna_hgvs:
        return render(template_file)

    # Invalid user input
    for hgvs in [dna_hgvs, rna_hgvs]:
        if error := validate_user_input(hgvs):
            # TODO make render work with
            return render(
                template_file,
                dna=dna_hgvs,
                rna=rna_hgvs,
                protein_domains=protein_domains,
                extended=extended,
                error=error,
            )

    # Analyze the transcript
    try:
        allele = Allele(dna=dna_hgvs, rna=rna_hgvs)

        if protein_domains:
            allele.lookup_protein_domains(allele.rna.d)

        results = allele.analyze(extended=extended)
    except Exception as e:
        error = {"summary": "Analysis failed", "details": str(e)}
        results = []

    # Get external links
    try:
        links = lookup_variant(dna_hgvs).url_dict()
    except Exception as e:
        links = dict()

    if error:
        return render(
            template_file,
            dna=dna_hgvs,
            rna=rna_hgvs,
            protein_domains=protein_domains,
            extended=extended,
            error=error,
        )
    else:
        return render(
            template_file,
            dna=dna_hgvs,
            rna=rna_hgvs,
            protein_domains=protein_domains,
            extended=extended,
            results=results,
            links=links,
            variant=rna_hgvs,
        )
