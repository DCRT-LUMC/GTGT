import gzip
import json
import logging
import os
import urllib.request
from abc import ABC, abstractmethod
from typing import Any, Protocol
from urllib.error import HTTPError

payload = dict[str, Any]
parameters = tuple[Any, ...]

logger = logging.getLogger(__name__)


class Provider(Protocol):
    """Protocol for the Provider class

    There is also an abstrac base class _Provider, which implements file system
    caching of API calls and can be inherited from
    """

    def get(self, parameters: parameters) -> payload: ...


class _Provider(ABC):
    def __init__(self) -> None:
        cache = os.environ.get("GTGT_CACHE")
        name = type(self).__name__

        self.cache: str | None
        # Put the cache for each Provider in a separate folder
        if cache:
            self.cache = f"{cache}/{name}"
        else:
            self.cache = None

        # Ensure the cache folder exists
        if self.cache:
            os.makedirs(self.cache, exist_ok=True)

    @abstractmethod
    def get(self, parameters: parameters) -> payload:
        pass

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

    def _fname(self, parameters: parameters) -> str:
        """Generate file name for the cache based on the parameters"""
        return f"{self.cache}/{'_'.join(str(x) for x in parameters)}.json.gz"

    def _get(self, url: str, parameters: parameters) -> payload:
        """Get the requested data, from the filename or the url"""
        # If the cache is not enabled
        if not self.cache:
            return self._fetch_url(url)

        # Filename for the cached payload
        fname = self._fname(parameters)

        js: payload = dict()
        # If the payload is already in the cache
        if os.path.exists(fname):
            logger.info(f"Reading payload from {fname}")
            with gzip.open(fname, "rt") as fin:
                js = json.load(fin)
        else:
            # If the payload is not in the cache
            js = self._fetch_url(url)
            with gzip.open(fname, "wt") as fout:
                print(json.dumps(js), file=fout)
        return js


class MyGene(_Provider):
    def get(self, parameters: parameters) -> payload:
        ensembl_gene_id, *rest = parameters
        url = f"https://mygene.info/v3/gene/{ensembl_gene_id}?fields=uniprot"

        return self._get(url, parameters)


class VariantValidator(_Provider):
    def get(self, parameters: parameters) -> payload:
        prefix = "https://rest.variantvalidator.org/VariantValidator"
        suffix = "mane_select?content-type=application/json"

        assembly, variant, *rest = parameters
        if variant.startswith("ENS"):
            url = f"{prefix}/variantvalidator_ensembl/{assembly}/{variant}/{suffix}"
        else:
            url = f"{prefix}/variantvalidator/{assembly}/{variant}/{suffix}"

        return self._get(url, parameters)


class UCSC(_Provider):
    def get(self, parameters: parameters) -> payload:
        genome, chrom, start, end, track = parameters
        url = ";".join(
            (
                f"https://api.genome.ucsc.edu/getData/track?genome={genome}",
                f"chrom={chrom}",
                f"start={start}",
                f"end={end}",
                f"track={track}",
            )
        )
        return self._get(url, parameters)


class Ensembl(_Provider):
    def get(self, parameters: parameters) -> payload:
        transcript, *rest = parameters
        url = f"http://rest.ensembl.org/lookup/id/{transcript}?content-type=application/json"
        return self._get(url, parameters)
