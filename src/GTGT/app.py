import uvicorn as uvicorn
from fastapi import FastAPI

from .variant_validator import lookup_variant
from .provider import Provider
from .bed import Bed
from .models import TranscriptModel, BedModel
from .wrappers import lookup_transcript

from typing import Dict

app = FastAPI()
provider = Provider()


@app.get("/links/{variant}")
async def get_links(variant: str) -> Dict[str, str]:
    return lookup_variant(provider, variant).url_dict()


@app.get("/transcript/{transcript_id}")
async def get_transcript(transcript_id: str) -> TranscriptModel:
    return lookup_transcript(provider, transcript_id)
