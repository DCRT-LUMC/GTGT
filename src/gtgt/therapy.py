import dataclasses
import logging
import os

# Only import Allele for type checking
from typing import TYPE_CHECKING, Any, Iterator, Mapping, Sequence, TypeVar

if TYPE_CHECKING:
    from gtgt.transcript import Allele

from gtgt.mutalyzer import (
    get_exons,
)
from gtgt.variant import Variant, combine_variants_deletion
from gtgt.vulexmap import VulExMap, lookup_vulexmap, vulexmap_description

logger = logging.getLogger(__name__)

# Initialize the VulExMap data fetcher once
V = VulExMap(path=os.environ.get("GTGT_VulExMap"))


@dataclasses.dataclass
class Therapy:
    """Internal representation of a generated Therapy"""

    name: str
    description: str
    rna_variants: list[Variant]
    dna_variants: list[Variant]


@dataclasses.dataclass
class TherapyReport:
    """Class to store genetic therapy information to present to the user"""

    name: str
    hgvsc: str
    description: str
    dna_variants: Sequence[Variant]
    rna_variants: Sequence[Variant]
    figure: str | None = None
    hgvsr: str | None = None
    hgvsp: str | None = None

    @classmethod
    def from_dict(cls, dict: Mapping[str, Any]) -> "TherapyReport":
        """Create a Therapy object from a dict representation of a Therapy"""
        dna_variants = [Variant.from_dict(x) for x in dict["dna_variants"]]
        rna_variants = [Variant.from_dict(x) for x in dict["rna_variants"]]
        return cls(
            name=dict["name"],
            hgvsc=dict["hgvsc"],
            description=dict["description"],
            dna_variants=dna_variants,
            rna_variants=rna_variants,
            figure=dict.get("figure"),
            hgvsr=dict.get("hgvsr"),
            hgvsp=dict.get("hgvsp"),
        )

    @classmethod
    def from_therapy(cls, therapy: Therapy, allele: "Allele") -> "TherapyReport":
        """Generate a TherapyReport from a Therapy and Allele"""
        return cls(
            name=therapy.name,
            hgvsc=allele.dna.hgvs(),
            description=therapy.description,
            dna_variants=allele.dna.variants,
            rna_variants=allele.rna.variants,
            hgvsr=allele.rna.hgvs(),
            hgvsp=allele.rna.protein_hgvs(),
        )


T = TypeVar("T")


def sliding_window(items: Sequence[T], size: int = 1) -> Sequence[Sequence[T]]:
    adj: list[Sequence[T]] = list()
    for i in range(len(items) - size + 1):
        adj.append([x for x in items[i : i + size]])
    return adj


def _exon_string(exon_numbers: Sequence[int]) -> str:
    """Format the exon names for a variable number of exons"""
    if len(exon_numbers) == 1:
        return f"exon {exon_numbers[0]}"
    elif len(exon_numbers) == 2:
        return f"exons {exon_numbers[0]} and {exon_numbers[1]}"
    else:
        t = ", ".join((str(x) for x in exon_numbers[:-1]))
        return f"exons {t} and {exon_numbers[-1]}"


def skip_adjacent_exons(allele: "Allele", number_to_skip: int = 1) -> Iterator[Therapy]:
    """Skipp all possible adjacent exons the specified Description"""

    skippable_exons = get_exons(allele.rna.d, in_transcript_order=True)[1:-1]
    variants = allele.rna.variants

    for i, exons in enumerate(sliding_window(skippable_exons, size=number_to_skip), 2):
        # Generate the string of exon numbers
        exons_description = _exon_string(range(i, i + number_to_skip))

        if allele.rna.d.is_inverted():
            exons = exons[::-1]

        # Start of the first exon to skip
        start = exons[0][0]
        # End of the last exon to skip
        end = exons[-1][-1]

        exon_skip = Variant(start, end)

        # Combine the existing variants with the exon skip
        try:
            combined = combine_variants_deletion(variants, exon_skip)
        except ValueError as e:
            if number_to_skip == 1:
                msg = f"Cannot skip exon {exons_description}: {e}"
            else:
                msg = f"Cannot skip exons {exons_description}: {e}"
            logger.warn(msg)
            continue
        logger.debug(f"Skip {exons_description}({exons=}): {exon_skip=} {combined=}")

        description = [
            f"The annotations based on the supplied variants, in combination with skipping {exons_description}.\n"
        ]
        # Lookup VulExMap data fore every exon
        for j, exon in enumerate(exons):
            vulexmap = lookup_vulexmap(allele.rna.d, exon, V)
            if vulexmap is not None:
                description.append(vulexmap_description(vulexmap, name=f"Exon {i+j}"))
        # If we added VulExMap information
        if len(description) > 1:
            description.append(
                "Please see https://vulexmap.compbio.sdu.dk/ for more details."
            )

        yield Therapy(
            name=f"Skip {exons_description}",
            description=" ".join(description),
            dna_variants=allele.dna.variants,
            rna_variants=list(combined),
        )


def generate_therapies(allele: "Allele") -> Iterator[Therapy]:
    """Wrapper around the different therapies"""
    # Skip a single exon
    yield from skip_adjacent_exons(allele, number_to_skip=1)
    # Skip two adjacent exons
    yield from skip_adjacent_exons(allele, number_to_skip=2)
