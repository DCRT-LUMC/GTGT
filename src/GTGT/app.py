import uvicorn as uvicorn
from fastapi import FastAPI, Body

from .variant_validator import lookup_variant
from .provider import Provider
from .models import BedModel, TranscriptModel
from .wrappers import lookup_transcript

from typing import Dict
from typing_extensions import Annotated

app = FastAPI()
provider = Provider()


@app.get("/links/{variant}")
async def get_links(variant: Annotated[str, "NM_000094.4:c.5299G>C"]) -> Dict[str, str]:
    """Lookup external references for the specified variant"""
    return lookup_variant(provider, variant).url_dict()


@app.get("/transcript/{transcript_id}")
async def get_transcript(transcript_id: str) -> TranscriptModel:
    """Lookup the specified transcript"""
    return lookup_transcript(provider, transcript_id)


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
) -> Dict[str, float]:
    """Compare two transcripts"""
    s = self.to_transcript()
    o = other.to_transcript()

    return s.compare(o)
