import json
import logging
import os
import urllib.request
from abc import ABC, abstractmethod
from typing import Any
from urllib.error import HTTPError

payload = dict[str, Any]

logger = logging.getLogger(__name__)


class Provider(ABC):
    def __init__(self) -> None:
        cache = os.environ.get("GTGT_CACHE")
        name = type(self).__name__

        # Put the cache for each Provider in a separate folder
        if cache:
            self.cache = f"{cache}/{name}"

        # Ensure the cache folder exists
        if self.cache:
            os.makedirs(self.cache, exist_ok=True)

    def __str__(self) -> str:
        return f"{type(self).__name__}(cache={self.cache})"

    def _fetch_url(self, url: str) -> payload:
        logger.info(f"Fetching {url=}")
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

    @abstractmethod
    def get(self, *args: Any) -> payload:
        pass

    def _get(self, url: str, fname: str) -> payload:
        """Get the requested data, from the filename or the url"""
        # If the cache is not enabled
        if not self.cache:
            return self._fetch_url(url)

        js: payload = dict()
        # If the payload is already in the cache
        if os.path.exists(fname):
            logger.info(f"Reading payload from {fname}")
            with open(fname) as fin:
                js = json.load(fin)
        else:
            # If the payload is not in the cache
            js = self._fetch_url(url)
            with open(fname, "wt") as fout:
                print(json.dumps(js), file=fout)
        return js


class MyGene(Provider):
    def get(self, ensembl_gene_id: str) -> payload:
        url = f"https://mygene.info/v3/gene/{ensembl_gene_id}?fields=uniprot"
        fname = f"{self.cache}/{ensembl_gene_id}.json"

        return self._get(url, fname)


class VariantValidator(Provider):
    def get(self, assembly: str, variant: str) -> payload:
        prefix = "https://rest.variantvalidator.org/VariantValidator"
        suffix = "mane_select?content-type=application/json"

        if variant.startswith("ENS"):
            url = f"{prefix}/variantvalidator_ensembl/{assembly}/{variant}/{suffix}"
        else:
            url = f"{prefix}/variantvalidator/{assembly}/{variant}/{suffix}"

        fname = f"{self.cache}/{assembly}_{variant}.json"

        return self._get(url, fname)


class UCSC(Provider):
    def get(self, genome: str, chrom: str, start: int, end: int, track: str) -> payload:
        url = ";".join(
            (
                f"https://api.genome.ucsc.edu/getData/track?genome={genome}",
                f"chrom={chrom}",
                f"start={start}",
                f"end={end}",
                f"track={track}",
            )
        )
        fname = f"{self.cache}/{genome}_{chrom}:{start}-{end}_{track}.json"
        return self._get(url, fname)


class Ensembl(Provider):
    def get(self, transcript: str) -> payload:
        url = f"http://rest.ensembl.org/lookup/id/{transcript}?content-type=application/json"
        fname = f"{self.cache}/{transcript}.json"
        return self._get(url, fname)
