from GTGT.variant_validator import lookup_variant
from .provider import Provider
from .wrappers import lookup_transcript
from flask import Flask, render_template

from typing import Optional

app = Flask(__name__)
provider = Provider()


@app.route("/")
@app.route("/<variant>")
def result(variant: Optional[str] = None) -> str:
    template_file = "index.html.j2"
    if not variant:
        return render_template(template_file)

    # Analyze the transcript
    transcript_id = variant.split(":")[0]
    transcript_model = lookup_transcript(provider, transcript_id)
    transcript = transcript_model.to_transcript()
    results = transcript.analyze(variant)

    # Get external links
    links = lookup_variant(provider, variant).url_dict()

    return render_template(template_file, results=results, links=links, variant=variant)
