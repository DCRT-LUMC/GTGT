import uvicorn as uvicorn
from fastapi import FastAPI

from .variant_validator import lookup_variant
from .provider import Provider

from typing import Dict

app = FastAPI()
provider = Provider()


@app.get("/links/{variant}")
async def get_links(variant: str) -> Dict[str, str]:
    return lookup_variant(provider, variant).url_dict()
