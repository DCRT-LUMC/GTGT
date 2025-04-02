from .variant_validator import lookup_variant
from .provider import Provider
from .wrappers import lookup_transcript

from flask import Flask, flash, render_template
from typing import Optional

import mutalyzer_hgvs_parser

hgvs_error = (
    mutalyzer_hgvs_parser.exceptions.UnexpectedCharacter,
    mutalyzer_hgvs_parser.exceptions.UnexpectedEnd,
)

app = Flask(__name__)
provider = Provider()


@app.route("/")
@app.route("/<variant>")
def result(variant: Optional[str] = None) -> str:
    template_file = "index.html.j2"

    # If no variant was specified
    if not variant:
        return render_template(template_file)

    # Test if the variant is valid HGVS
    try:
        mutalyzer_hgvs_parser.to_model(variant)
    except hgvs_error as e:
        flash(str(e))
        print("ERROR IS")
        print(e)
        return render_template(template_file)

    # Analyze the transcript
    transcript_id = variant.split(":")[0]
    transcript_model = lookup_transcript(provider, transcript_id)
    transcript = transcript_model.to_transcript()
    results = transcript.analyze(variant)

    # Get external links
    links = lookup_variant(provider, variant).url_dict()

    return render_template(template_file, results=results, links=links, variant=variant)
