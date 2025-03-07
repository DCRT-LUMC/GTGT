from GTGT.variant_validator import lookup_variant
from .provider import Provider
from .wrappers import lookup_transcript
from flask import Flask, render_template

app = Flask(__name__)
provider = Provider()

@app.route("/")
def home():
    return render_template("home.html")

@app.route("/<variant>")
def result(variant):
    if not variant:
        return render_template("report.html")

    # Analyze the transcript
    transcript_id = variant.split(":")[0]
    transcript_model = lookup_transcript(provider, transcript_id)
    transcript = transcript_model.to_transcript()
    results = transcript.analyze(variant)

    # Get external links
    links = lookup_variant(provider, variant).url_dict()

    print(results)
    return render_template("report.html", results=results, links=links, variant=variant)
