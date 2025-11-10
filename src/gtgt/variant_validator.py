from typing import Any, Mapping, Sequence, cast

from pydantic import BaseModel

from .provider import Provider

Payload = Mapping[str, Any]


class Links(BaseModel):
    omim_ids: Sequence[str]
    gene_symbol: str
    ensembl_gene_id: str
    uniprot: str
    decipher: str
    genomic_variant: str
    hgnc: str
    ucsc: str

    databases: Sequence[str] = [
        "omim",
        "lovd",
        "gtex",
        "uniprot",
        "decipher",
        "clinvar",
        "hgnc",
        "ucsc",
        "gnomad",
        "stringdb",
    ]

    def url(self, field: str) -> str | Sequence[str]:
        if field == "omim":
            urls = list()
            for id in self.omim_ids:
                url = f"https://www.omim.org/entry/{id}"
                urls.append(url)
            return urls
        elif field == "lovd":
            return f"https://databases.lovd.nl/shared/genes/{self.gene_symbol}"
        elif field == "gtex":
            return f"https://gtexportal.org/home/gene/{self.ensembl_gene_id}"
        elif field == "uniprot":
            return f"https://www.uniprot.org/uniprotkb/{self.uniprot}/entry"
        elif field == "decipher":
            return f"https://www.deciphergenomics.org/sequence-variant/{self.decipher}"
        elif field == "clinvar":
            return f"https://www.ncbi.nlm.nih.gov/clinvar/?term={self.genomic_variant}"
        elif field == "hgnc":
            return f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{self.hgnc}"
        elif field == "ucsc":
            return f"https://genome.cse.ucsc.edu/cgi-bin/hgGene?hgg_gene={self.ucsc}"
        elif field == "gnomad":
            return f"https://gnomad.broadinstitute.org/variant/{self.decipher}?dataset=gnomad_r4"
        elif field == "stringdb":
            return f"https://string-db.org/cgi/network?identifiers={self.gene_symbol}"
        else:
            raise NotImplementedError(f"Unknown field: '{field}'")

    def url_dict(self) -> Mapping[str, str]:
        """Create a flat dict with urls to all databases"""
        d = dict()

        for field in self.databases:
            # omim can contain a list of IDs
            if field == "omim":
                for i, url in enumerate(self.url(field), 1):
                    d[f"{field}_{i}"] = url
            else:
                d[field] = cast(str, self.url(field))

        return d


def guess_refseq_ensembl(variant: str) -> str:
    """Guess if a transcript is RefSeq or Ensembl"""
    if variant.startswith("ENS"):
        return "ensembl"
    else:
        return "refseq"


def lookup_variant(provider: Provider, variant: str, assembly: str = "hg38") -> Links:
    if guess_refseq_ensembl(variant) == "ensembl":
        url = f"https://rest.variantvalidator.org/VariantValidator/variantvalidator_ensembl/{assembly}/{variant}/mane_select?content-type=application/json"
    else:
        url = f"https://rest.variantvalidator.org/VariantValidator/variantvalidator/{assembly}/{variant}/mane_select?content-type=application/json"

    payload = provider.get(url)

    d = parse_payload(payload, variant, assembly)
    d["uniprot"] = lookup_uniprot(provider, d["ensembl_gene_id"])

    return Links(**d)


def lookup_uniprot(provider: Provider, ensembl_gene_id: str) -> str:
    url = f"https://mygene.info/v3/gene/{ensembl_gene_id}?fields=uniprot"
    uniprot_id: str = provider.get(url)["uniprot"]["Swiss-Prot"]
    return uniprot_id


def extract_variant(payload: Payload, variant: str) -> Payload:
    """Extract the variant section from the payload"""
    for value in payload.values():
        # Skip flag field
        if not isinstance(value, dict):
            continue
        submitted_variant = value.get("submitted_variant")
        if submitted_variant == variant:
            var_payload: Payload = value
            return var_payload
    else:
        msg = f"Unable to parse VariantValidator output, '{variant}' not found."
        raise ValueError(msg)


def parse_payload(payload: Payload, variant: str, assembly: str) -> dict[str, Any]:
    # Check the flag to see if the reply is valid
    flag = payload["flag"]
    if flag == "warning":
        w = payload.get("validation_warning_1", dict())
        errors = "\n".join(w.get("validation_warnings", []))
        raise ValueError(errors)
    if flag != "gene_variant":
        msg = f"Unknown VariantValidator flag: {flag}"
        raise NotImplementedError(msg)

    var = extract_variant(payload, variant)
    d = {
        "omim_ids": var["gene_ids"]["omim_id"],
        "gene_symbol": var["gene_symbol"],
        "ensembl_gene_id": var["gene_ids"]["ensembl_gene_id"],
        "hgnc": var["annotations"]["db_xref"]["hgnc"],
        "ucsc": var["gene_ids"]["ucsc_id"],
        "genomic_variant": var["primary_assembly_loci"][assembly][
            "hgvs_genomic_description"
        ],
    }
    # Get the 'VCF' notation on the specified assembly
    vcf = var["primary_assembly_loci"][assembly]["vcf"]
    # Remove the 'chr' prefix
    vcf["chr"] = vcf["chr"][3:]
    # Decypher uses {chrom}-{pos}-{ref}-{alt} as variant ID
    d["decipher"] = "-".join(vcf[field] for field in ["chr", "pos", "ref", "alt"])
    return d
