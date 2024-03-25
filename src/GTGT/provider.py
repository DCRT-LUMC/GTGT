from typing import Any, Dict, Optional
import json
import urllib.request
from urllib.error import HTTPError
from urllib.parse import urlparse

import logging

payload = Dict[str, Any]

logger = logging.getLogger(__name__)


class Provider:
    def __init__(self, cache_dir: Optional[str] = None):
        self.cache_dir = cache_dir

    def get(self, url: str) -> payload:
        logger.debug(f"Fetching {url=}")
        try:
            response = urllib.request.urlopen(url)
        except HTTPError as e:
            raise RuntimeError(e)
        else:
            js: Dict[str, Any] = json.loads(response.read())

        return js

    def url_to_filename(self, url: str) -> str:
        fname = ""
        parsed = urlparse(url)
        fname += parsed.netloc
        if len(parsed.path) > 1:
            fname += parsed.path.replace("/", "_")

        if parsed.query:
            fname += "_" + parsed.query.replace("/", "_")

        return fname
