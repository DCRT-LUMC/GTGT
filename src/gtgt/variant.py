from typing import Any, Mapping, Sequence

from mutalyzer.description import Description
from mutalyzer.description_model import get_reference_id
from schema import And, Optional, Or, Schema


def get_offset(d: Description) -> int:
    ref_id = get_reference_id(d.corrected_model)
    offset: int = (
        d.references.get(ref_id, {})
        .get("annotations", {})
        .get("qualifiers", {})
        .get("location_offset", 0)
    )
    return offset


class Variant:
    """Class to store delins variants"""

    # fmt: off
    # Schema for the location specification of the indel model
    location_schema = Schema(
        {
            "type": "range",
            "start": {
                "type": "point",
                "position": int,
            },
            "end": {
                "type": "point",
                "position": int,
            },
        }
    )

    # Schema for the inserted/deleted entries of the indel model
    inserted_deleted_schema = Schema(
        And( # Inserted must be 0 or 1 items
            lambda n: len(n) <= 1,
            [
                {
                    "sequence": Or(str, []),
                    "source": "description",
                    Optional("inverted") : True
                },
            ],
        ),
    )

    # Full schema for the indel model
    schema = Schema(
        {
            "type": "deletion_insertion",
            "source": "reference",
            "location": location_schema,
            Optional("inserted"): inserted_deleted_schema,
            Optional("deleted"): inserted_deleted_schema,
        }
    )
    # fmt: on

    def __init__(
        self,
        start: int,
        end: int,
        inserted: str = "",
        deleted: str = "",
    ):
        if start > end:
            raise ValueError(f"End ({end}) must be after start ({start})")
        self.start = start  # zero based
        self.end = end  # exclusive
        self.inserted = inserted

        if len(deleted) > 1:
            raise ValueError("deleted sequence is only defined for SNPS, not indels")
        self.deleted = deleted

    def __str__(self) -> str:
        return self.__repr__()

    def __repr__(self) -> str:
        start = self.start
        end = self.end
        inserted = self.inserted
        deleted = self.deleted
        return f"Variant({start=}, {end=}, inserted={inserted}, deleted={deleted})"

    def __hash__(self) -> int:
        return hash(self.__repr__())

    def before(self, other: "Variant") -> bool:
        return self.end <= other.start

    def after(self, other: "Variant") -> bool:
        return self.start >= other.end

    def inside(self, other: "Variant") -> bool:
        return self.start >= other.start and self.end <= other.end

    def overlap(self, other: "Variant") -> bool:
        self_ends_in_other = self.end > other.start and self.end <= other.end
        self_starts_in_other = self.start >= other.start and self.start < other.end

        return any(
            [
                self_starts_in_other,
                self_ends_in_other,
                self.inside(other),
                other.inside(self),
            ]
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Variant):
            raise NotImplementedError

        return (
            self.start == other.start
            and self.end == other.end
            and self.inserted == other.inserted
            and self.deleted == other.deleted
        )

    def __lt__(self, other: "Variant") -> bool:
        if not isinstance(other, Variant):
            raise NotImplementedError
        if self.overlap(other):
            msg = f"Overlapping variants '{self}' and '{other}' cannot be sorted"
            raise ValueError(msg)

        return self.start < other.start

    @staticmethod
    def _validate_schema(model: Mapping[str, Any]) -> None:
        """Validate the structure of the mutalyzer delins model

        This can be very complex, and we only support the most common cases.
        """

        Variant.schema.validate(model)

    @staticmethod
    def _model_is_repeat(model: Mapping[str, Any]) -> bool:
        """Determine if model is a repeat"""
        # Determine if the model is a repeat
        repeat_schema = Schema(
            And(  # An empty list would be a deletion, not a repeat
                lambda n: len(n) == 1,
                [
                    {
                        "source": "description",
                        "sequence": str,
                        "repeat_number": {"type": "point", "value": int},
                        Optional("inverted"): True,
                    }
                ],
            )
        )

        inserted = model.get("inserted")
        is_repeat: bool = repeat_schema.is_valid(inserted)

        return is_repeat

    @staticmethod
    def _model_repeat_to_delins(model: Mapping[str, Any]) -> dict[str, Any]:
        """Convert a repeat model to a delins model"""
        new_model = {k: v for k, v in model.items()}

        # Determine the sequence and repeats
        old_sequence = model["inserted"][0]["sequence"]
        repeats = model["inserted"][0]["repeat_number"].get("value", 1)

        # Expand the new sequence
        new_sequence = old_sequence * repeats

        inserted = {"sequence": new_sequence, "source": "description"}

        # If the model is defined on the reverse strand
        inverted = model["inserted"][0].get("inverted", False)
        if inverted:
            inserted["inverted"] = True

        new_model["inserted"] = [inserted]
        return new_model

    @staticmethod
    def _model_is_duplication(model: Mapping[str, Any]) -> bool:
        """Determine if model is a duplication"""
        # Determine if the model is a repeat
        location_schema = Schema(
            {
                "type": "range",
                "start": {
                    "type": "point",
                    "position": int,
                },
                "end": {
                    "type": "point",
                    "position": int,
                },
            }
        )
        duplication_schema = Schema(
            And(  # An empty list would be a deletion, not a repeat
                lambda n: len(n) == 1,
                [
                    {
                        "source": "reference",
                        "location": location_schema,
                        "repeat_number": {"type": "point", "value": int},
                    }
                ],
            )
        )

        inserted = model.get("inserted")
        is_duplication: bool = duplication_schema.is_valid(inserted)

        return is_duplication

    @staticmethod
    def _model_duplication_to_delins(
        model: Mapping[str, Any], sequence: str
    ) -> dict[str, Any]:
        """Convert a duplication model to a delins model"""
        if not sequence:
            raise ValueError("Variant: specify sequence to handle duplications")
        new_model = {k: v for k, v in model.items()}

        # Determine the start and end on the sequence
        inserted = model["inserted"][0]
        start = inserted["location"]["start"]["position"]
        end = inserted["location"]["end"]["position"]
        repeats = inserted["repeat_number"].get("value", 1)

        # Expand the new sequence
        new_sequence = sequence[start:end]

        new_model["inserted"] = [
            {
                "sequence": new_sequence,
                "repeat_number": {"type": "point", "value": repeats},
                "source": "description",
            }
        ]
        return new_model

    @staticmethod
    def _model_is_inversion(model: Mapping[str, Any]) -> bool:
        """Determine if model is an inversion"""
        # Determine if the model is a repeat
        location_schema = Schema(
            {
                "type": "range",
                "start": {
                    "type": "point",
                    "position": int,
                },
                "end": {
                    "type": "point",
                    "position": int,
                },
            }
        )
        duplication_schema = Schema(
            And(  # An empty list would be a deletion, not a repeat
                lambda n: len(n) == 1,
                [
                    {
                        "source": "reference",
                        "location": location_schema,
                        "inverted": True,
                    }
                ],
            )
        )

        inserted = model.get("inserted")
        is_duplication: bool = duplication_schema.is_valid(inserted)

        return is_duplication

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        return "".join((complement[nt] for nt in reversed(seq)))

    @staticmethod
    def _model_inversion_to_delins(
        model: Mapping[str, Any], sequence: str
    ) -> dict[str, Any]:
        """Convert an inversion model to a delins model"""
        if not sequence:
            raise ValueError("Variant: specify sequence to handle inversions")
        new_model = {k: v for k, v in model.items()}

        # Determine the start and end on the sequence
        inserted = model["inserted"][0]
        start = inserted["location"]["start"]["position"]
        end = inserted["location"]["end"]["position"]

        # Expand the new sequence
        new_sequence = Variant._reverse_complement(sequence[start:end])

        new_model["inserted"] = [
            {
                "sequence": new_sequence,
                "source": "description",
            }
        ]
        return new_model

    @classmethod
    def from_model(cls, model: Mapping[str, Any], sequence: str = "") -> "Variant":
        if Variant._model_is_inversion(model):
            model = Variant._model_inversion_to_delins(model, sequence)
        if Variant._model_is_duplication(model):
            model = Variant._model_duplication_to_delins(model, sequence)
        # Determine if the model is a repeat
        if Variant._model_is_repeat(model):
            model = Variant._model_repeat_to_delins(model)

        # Validate the delins model
        cls._validate_schema(model)

        start = model["location"]["start"]["position"]
        end = model["location"]["end"]["position"]

        # Store if the inserted and deleted sequences were inverted
        ins_inverted: bool | None = None
        del_inverted: bool | None = None

        inserted = model.get("inserted", [])

        # Sanity check, complex variants should have been normalized into a
        # delins representation
        if len(inserted) == 0:
            inserted = ""
        elif len(inserted) > 1:
            raise NotImplementedError(
                "Multiple records in 'inserted' are not supported"
            )
        else:  # inserted contains 1 item, as is expected
            if "sequence" not in inserted[0]:
                raise NotImplementedError(
                    "Missing sequence in 'inserted' is not supported"
                )
            if "repeat_number" in inserted[0]:
                raise NotImplementedError("Repeats in 'inserted' are not supported")

            ins_inverted = inserted[0].get("inverted", False)
            inserted = inserted[0]["sequence"]

        deleted = model.get("deleted")
        if deleted is not None:
            try:
                del_inverted = deleted[0].get("inverted", False)
                deleted = deleted[0]["sequence"]
            except KeyError:
                raise NotImplementedError("Complex Variant not supported")

        # Make sure inserted and deleted are both or neither inversed
        if ins_inverted is not None and del_inverted is not None:
            if ins_inverted != del_inverted:
                msg = "strand difference between inserted and deleted sequences are not supported"
                raise NotImplementedError(msg)

        # If inserted or deleted are inverted, we need to reverse complement
        # the inserted/deleted nucleotides
        if inserted and ins_inverted:
            inserted = cls._reverse_complement(inserted)
        if deleted and del_inverted:
            deleted = cls._reverse_complement(deleted)

        return Variant(
            start=start,
            end=end,
            inserted=inserted if inserted else "",
            deleted=deleted if deleted else "",
        )

    @classmethod
    def from_dict(cls, dict: Mapping[str, Any]) -> "Variant":
        """Create a Variant object from a dict representation of a Variant"""
        return cls(**dict)

    def to_model(self) -> Mapping[str, Any]:
        """Convert Variant to mutalyzer delins model"""

        # Specification of the location
        # fmt: off
        location = {
            "type": "range",
            "start": {
                "type": "point",
                "position": self.start
            },
            "end": {
                "type": "point",
                "position": self.end
            }
        }

        # Specification of the inserted sequence
        inserted_obj: dict[str, Any] = {
            "sequence": self.inserted,
            "source": "description"
        }

        if self.inserted:
            inserted = [inserted_obj]
        else:
            inserted = []
        # fmt: on

        model = {
            "location": location,
            "type": "deletion_insertion",
            "source": "reference",
            "inserted": inserted,
        }

        deletion_object: dict[str, Any] = {
            "sequence": self.deleted,
            "source": "description",
        }
        if self.deleted:
            model["deleted"] = [deletion_object]

        return model

    def genomic_coordinates(self, d: Description) -> tuple[int, int]:
        """Return genomic coordinates for Variant"""
        offset = get_offset(d)

        if offset is None:
            raise RuntimeError("Missing ensembl offset")

        return self.start + offset, self.end + offset


def combine_variants_deletion(
    variants: Sequence[Variant], deletion: Variant
) -> Sequence[Variant]:
    """Combine variants and a deletion, any variants that are contained in
    the deletion are discarded

    The resulting list of variants is sorted
    """
    if deletion.inserted:
        raise ValueError(f"{Variant} is not a pure deletion")

    # Ensure the variants are sorted, and do not overlap
    sorted_variants = sorted(variants)

    combined = list()
    for i in range(len(sorted_variants)):
        variant = sorted_variants[i]
        if variant.before(deletion):
            combined.append(variant)
        elif variant.inside(deletion):
            # Discard the current variant
            continue
        elif variant.after(deletion):
            combined.append(deletion)
            combined += sorted_variants[i:]
            break
        else:
            msg = f"Deletion '{deletion}' partially overlaps '{variant}"
            raise ValueError(msg)
    else:
        combined.append(deletion)

    return combined
