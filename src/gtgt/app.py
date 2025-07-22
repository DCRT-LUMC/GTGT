from typing import Mapping, Sequence

import uvicorn as uvicorn
from fastapi import Body, FastAPI, HTTPException
from fastapi.responses import RedirectResponse
from typing_extensions import Annotated

from .models import BedModel, TranscriptId, TranscriptModel
from .mutalyzer import HGVS
from .provider import Provider
from .transcript import Comparison, Result
from .variant_validator import lookup_variant
from .wrappers import lookup_transcript

app = FastAPI()
provider = Provider()


@app.get("/")
async def redirect() -> RedirectResponse:
    response = RedirectResponse(url="/docs")
    return response


@app.post("/links")
async def get_links(variant: HGVS) -> Mapping[str, str]:
    """Lookup external references for the specified variant"""
    try:
        reply = lookup_variant(provider, variant.description).url_dict()
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
    return reply


@app.post("/transcript")
async def get_transcript(transcript: TranscriptId) -> TranscriptModel:
    """Lookup the specified transcript"""
    try:
        ts = lookup_transcript(provider, transcript.id)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
    return ts


@app.post("/transcript/exonskip")
async def exon_skip(transcript: TranscriptModel, region: BedModel) -> TranscriptModel:
    """Skip exons that overlap the specified region"""
    ts = transcript.to_transcript()
    skip_region = region.to_bed()
    ts.exon_skip(skip_region)
    return TranscriptModel.from_transcript(ts)


@app.post("/transcript/compare")
async def compare(
    self: Annotated[
        TranscriptModel,
        Body(
            examples=[
                {
                    "exons": {
                        "chrom": "chr1",
                        "blocks": [[0, 10], [50, 60], [70, 100]],
                        "name": "exons",
                        "score": 0,
                        "strand": ".",
                    },
                    "cds": {
                        "chrom": "chr1",
                        "blocks": [[40, 72]],
                        "name": "cds",
                        "score": 0,
                        "strand": ".",
                    },
                }
            ]
        ),
    ],
    other: TranscriptModel,
) -> Sequence[Comparison]:
    """Compare two transcripts"""
    s = self.to_transcript()
    o = other.to_transcript()

    return s.compare(o)


@app.post("/hgvs/analyze", response_model=None)
async def analyze(hgvs: HGVS) -> Sequence[Result]:
    """Analyze all possible exons skips for the spcified HGVS variant"""
    transcript_id = hgvs.description.split(":")[0]
    transcript_model = lookup_transcript(provider, transcript_id)
    transcript = transcript_model.to_transcript()
    return transcript.analyze(hgvs.description)
