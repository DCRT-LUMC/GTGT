from collections import defaultdict
from typing import Any, Mapping, Sequence, cast

from pydantic import BaseModel

from .provider import MyGene, Provider, VariantValidator

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

    databases: Mapping[str, str] = {
        "omim": "Gene",
        "lovd": "Variant",
        "gtex": "Gene",
        "uniprot": "Protein",
        "decipher": "Variant",
        "clinvar": "Variant",
        "hgnc": "Gene",
        "ucsc": "Gene",
        "gnomad": "Variant",
        "stringdb": "Protein",
    }

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

    def description(self, field: str) -> str:
        d = {
            "lovd": (
                "LOVD is a database of genetic variants organized by gene. Annotations "
                "include variant description on the transcript and genome level, as well as "
                "the clinical classification of each variant."
            ),
            "gtex": (
                "GTEx contains the expression levels of a gene across different "
                "tissues derived from bulk sequencing data as well as from single cell experiments."
            ),
            "uniprot": (
                "The Universal Protein Resource (UniProt) is a comprehensive "
                "resource for protein sequence and annotation data. "
                "The mission of UniProt is to provide the scientific "
                "community with a comprehensive, high-quality and freely "
                "accessible resource of protein sequence and functional information."
            ),
            "decipher": (
                "DECIPHER (DatabasE of genomiC varIation and Phenotype in "
                "Humans using Ensembl Resources) is an interactive web-based "
                "database which incorporates a suite of tools designed to aid "
                "the interpretation of genomic variants. "
                "DECIPHER enhances clinical diagnosis by retrieving "
                "information from a variety of bioinformatics resources "
                "relevant to the variant found in the patient. The patient's "
                "variant is displayed in the context of both normal variation "
                "and pathogenic variation reported at that locus thereby "
                "facilitating interpretation."
            ),
            "clinvar": (
                "ClinVar is a freely accessible, public archive of reports of "
                "human variations classified for diseases and drug responses, "
                "with supporting evidence. ClinVar thus facilitates access to "
                "and communication about the relationships asserted between "
                "human variation and observed conditions, and the history of "
                "those assertions"
            ),
            "hgnc": (
                "The HGNC is responsible for approving unique symbols and "
                "names for human loci, including protein coding genes, ncRNA "
                "genes and pseudogenes, to allow unambiguous scientific communication"
            ),
            "ucsc": (
                "The Genome Browser stacks annotation tracks beneath genome "
                "coordinate positions, allowing rapid visual correlation of "
                "different types of information. The user can look at a whole "
                "chromosome to get a feel for gene density, open a specific "
                "cytogenetic band to see a positionally mapped disease gene "
                "candidate, or zoom in to a particular gene to view its "
                "spliced ESTs and possible alternative splicing."
            ),
            "gnomad": (
                "The Genome Aggregation Database (gnomAD™), originally "
                "launched in 2014 as the Exome Aggregation Consortium (ExAC), "
                "is the result of a coalition of investigators willing to share "
                "aggregate human exome and genome sequencing data from a "
                "variety of large-scale sequencing projects, and make summary "
                "data available for the wider scientific community."
            ),
            "stringdb": (
                "STRING represents relationships between proteins as a network "
                "(graph). In the network visualization, proteins are shown as "
                "nodes (bubbles) and associations between them are shown as "
                "edges (lines). Each node represents a protein encoded by a "
                "single gene locus."
            ),
            "omim": (
                "Online Mendelian Inheritance in Man®. OMIM is a comprehensive, "
                "authoritative compendium of human genes and genetic phenotypes "
                "that is freely available and updated daily. The full-text, "
                "referenced overviews in OMIM contain information on all known "
                "mendelian disorders and over 16,000 genes"
            ),
        }
        if field.startswith("omim"):
            field = "omim"
        return d.get(field, f"Unknown resource ({field})")

    def url_dict(self) -> Mapping[str, Any]:
        """Create a nested dict with urls to all databases"""
        d = defaultdict(list)

        for field, category in self.databases.items():
            # omim can contain a list of IDs
            if field == "omim":
                for i, url in enumerate(self.url(field), 1):
                    entry = {
                        "name": f"{field}_{i}",
                        "url": url,
                        "description": self.description(field),
                    }
                    d[category].append(entry)
            else:
                entry = {
                    "name": field,
                    "url": cast(str, self.url(field)),
                    "description": self.description(field),
                }
                d[category].append(entry)

        return d


def lookup_variant(variant: str, assembly: str = "hg38") -> Links:
    provider: Provider = VariantValidator()
    payload = provider.get((assembly, variant))

    d = parse_payload(payload, variant, assembly)
    d["uniprot"] = lookup_uniprot(d["ensembl_gene_id"])

    return Links(**d)


def lookup_uniprot(ensembl_gene_id: str) -> str:
    provider: Provider = MyGene()

    payload = provider.get((ensembl_gene_id,))
    uniprot_id: str = payload["uniprot"]["Swiss-Prot"]
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
