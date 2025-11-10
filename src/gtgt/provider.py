import json
import logging
import os
import urllib.request
from typing import Any
from urllib.error import HTTPError
from urllib.parse import urlparse

payload = dict[str, Any]

logger = logging.getLogger(__name__)


class Provider:
    def __init__(self, cache_dir: str | None = None):
        self.cache_dir = cache_dir

    def get(self, url: str) -> payload:
        if self.cache_dir is None:
            return self._fetch_url(url)
        else:
            return self._fetch_cache(url)

    def _fetch_url(self, url: str) -> payload:
        logger.debug(f"Fetching {url=}")
        try:
            response = urllib.request.urlopen(url)
        except HTTPError as e:
            raise RuntimeError(str(e))

        data = response.read()

        try:
            js: payload = json.loads(data)
        except Exception as e:
            logger.error(data)
            raise e

        return js

    def _fetch_cache(self, url: str) -> payload:
        """Fetch cache based on the provided url

        If no cache exists, use self._fetch_url to get and store the data
        """
        fname = f"{self.cache_dir}/{self.url_to_filename(url)}"
        data: payload = dict()

        # Try to open the cache file
        if os.path.exists(fname):
            with open(fname) as fin:
                data = json.load(fin)
                logger.debug(f"Retrieved cache for {url=}")
        else:
            data = self._fetch_url(url)

            # Ensure the folder exists
            os.makedirs(self.cache_dir, exist_ok=True)  # type: ignore

            with open(fname, "w") as fout:
                logger.debug(f"Writing cache for {url=}")
                print(json.dumps(data, indent=True), file=fout)
        return data

    def url_to_filename(self, url: str) -> str:
        fname = ""
        parsed = urlparse(url)
        fname += parsed.netloc
        if len(parsed.path) > 1:
            fname += parsed.path.replace("/", "_")

        if parsed.query:
            fname += "_" + parsed.query.replace("/", "_")

        return fname
