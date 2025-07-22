import pytest

from gtgt.provider import Provider


def test_default_provider() -> None:
    provider = Provider()
    assert provider.cache_dir is None


URLS = [
    ("http://rest.ensembl.org/", "rest.ensembl.org"),
    (
        "http://rest.ensembl.org/lookup/ENST00000241453",
        "rest.ensembl.org_lookup_ENST00000241453",
    ),
    (
        "http://rest.ensembl.org/lookup/ENST00000241453?content-type=application/json",
        "rest.ensembl.org_lookup_ENST00000241453_content-type=application_json",
    ),
    (
        "https://api.ucsc.edu/getData/track?genome=hg38;chrom=chr13;track=knownGene;start=28003274;end=28100576",
        "api.ucsc.edu_getData_track_genome=hg38;chrom=chr13;track=knownGene;start=28003274;end=28100576",
    ),
]


@pytest.mark.parametrize("url, expected", URLS)
def test_url_to_filename(url: str, expected: str) -> None:
    P = Provider("cache")
    assert P.url_to_filename(url) == expected
