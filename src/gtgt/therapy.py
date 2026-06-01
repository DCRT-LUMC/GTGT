import dataclasses
import logging
import os
from typing import Any, Mapping, Sequence, TypeVar

from mutalyzer.description import Description

from gtgt.mutalyzer import (
    get_exons,
    protein_prediction,
    sequence_from_description,
    to_cdot_hgvs,
)
from gtgt.variant import Variant, combine_variants_deletion
from gtgt.vulexmap import VulExMap, lookup_vulexmap, vulexmap_description

logger = logging.getLogger(__name__)

# Initialize the VulExMap data fetcher once
V = VulExMap(path=os.environ.get("GTGT_VulExMap"))


@dataclasses.dataclass
class Therapy:
    """Class to store genetic therapies"""

    name: str
    hgvsc: str
    description: str
    variants: Sequence[Variant]
    figure: str | None = None
    hgvsr: str | None = None
    hgvsp: str | None = None

    @classmethod
    def from_dict(cls, dict: Mapping[str, Any]) -> "Therapy":
        """Create a Therapy object from a dict representation of a Therapy"""
        v = [Variant.from_dict(x) for x in dict["variants"]]
        return cls(
            name=dict["name"],
            hgvsc=dict["hgvsc"],
            description=dict["description"],
            variants=v,
            figure=dict.get("figure"),
            hgvsr=dict.get("hgvsr"),
            hgvsp=dict.get("hgvsp"),
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


def skip_adjacent_exons(d: Description, number_to_skip: int = 1) -> Sequence[Therapy]:
    """Skipp all possible adjacent exons the specified Description"""
    exon_skips: list[Therapy] = list()

    skippable_exons = get_exons(d, in_transcript_order=True)[1:-1]
    sequence = sequence_from_description(d)
    variants = [
        Variant.from_model(v, sequence=sequence) for v in d.delins_model["variants"]
    ]
    logger.debug(f"Input variants: {variants}")

    for i, exons in enumerate(sliding_window(skippable_exons, size=number_to_skip), 2):
        # Generate the string of exon numbers
        exons_description = _exon_string(range(i, i + number_to_skip))

        if d.is_inverted():
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
            vulexmap = lookup_vulexmap(d, exon, V)
            if vulexmap is not None:
                description.append(vulexmap_description(vulexmap, name=f"Exon {i+j}"))
        # If we added VulExMap information
        if len(description) > 1:
            description.append(
                "Please see https://vulexmap.compbio.sdu.dk/ for more details."
            )

        # Convert to c. notation (user facing)
        name = f"Skip {exons_description}"
        selector = d.get_selector_id()
        cdot_variants = to_cdot_hgvs(d, combined)
        hgvsc = f"{selector}:c.{cdot_variants}"
        hgvsr = f"{selector}:r.{cdot_variants}"
        description = description
        t = Therapy(
            name=name,
            hgvsc=hgvsc,
            hgvsr=hgvsr,
            hgvsp=protein_prediction(d, combined)[0],
            description=" ".join(description),
            variants=combined,
        )
        exon_skips.append(t)

    return exon_skips


def generate_therapies(d: Description) -> Sequence[Therapy]:
    """Wrapper around the different therapies"""
    therapies: list[Therapy] = list()
    # Skip a single exon
    therapies += skip_adjacent_exons(d, number_to_skip=1)
    # Skip two adjacent exons
    therapies += skip_adjacent_exons(d, number_to_skip=2)
    return therapies
