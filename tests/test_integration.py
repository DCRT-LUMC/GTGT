from typing import Sequence

import pytest

from gtgt import Transcript
from gtgt.mutalyzer import init_description, sequence_from_description
from gtgt.transcript import Comparison
from gtgt.variant import Variant


def compare_to_wildtype(hgvs: str) -> Sequence[Comparison]:
    """
    Helper function to build a transcript from the supplied HGVS, and compare it to the wildtype transcript
    """
    # First, create the wildtype transcript
    *rest, variant = hgvs.split(".")
    wt_hgvs = ".".join(rest) + ".="

    wt = Transcript.from_description(init_description(wt_hgvs))

    # Next, we create and mutate the transcript of interest
    d = init_description(hgvs)
    t = Transcript.from_description(d)

    sequence = sequence_from_description(d)
    input_variants = [
        Variant.from_model(delins, sequence=sequence)
        for delins in d.delins_model["variants"]
    ]
    t.mutate(d, input_variants)

    return t.compare(wt)


class TestDifferentTranscripts:
    """Class to test the same mutation on different transcripts
    (with different internal coordinate systems)
    """

    coordinate = ("coordinate", ["c", "r"])

    SDHD_transcript = (
        "transcript",
        ["ENST00000375549.8", "NC_000011.10(NM_003002.4)", "NM_003002.4"],
    )
    WT1_transcript = (
        "transcript",
        ["ENST00000452863.10", "NC_000011.test(NM_024426.6)", "NM_024426.6"],
    )

    # NOTE: These variants are special, they are all deletions of one or more full amino acids
    SDHD_variants = [
        # Deletion in exon 1
        "31_33del",
        # Deletion across exon 1 and 2
        "49_57del",
        # Deletion in exon 2
        "100_132del",
        # Deletion in the last exon
        "451_453del",
        # Multiple variants
        "[31_33del;451_453del]",
    ]

    # NOTE: These variants are special, they are all deletions of one or more full amino acids
    WT1_variants = [
        # Deletion in exon 1
        "4_6del",
        "7_9del",
        # Deletion crosses exon 1 and exon2 boundary
        "658_663del",
        # Deletion in exon 2 (not normalized, which whould shift the variant 1bp downstream)
        "700_720del",
        # Deletion in the last exon (not normalized)
        "1252_1260del",
        # Multiple variants
        "[7_9del;658_663del]",
    ]

    @pytest.mark.slow
    @pytest.mark.parametrize(*coordinate)
    @pytest.mark.parametrize("variant", SDHD_variants)
    def test_different_transcripts_SDHD(self, coordinate: str, variant: str) -> None:
        """Test that ENST, NC(NM) and NM variants are handled the same

        Because the internal coordinate systems differ we cannot compare these
        directly but we have to compare them via the comparison with the wildtype
        version of the corresponding transcript

        Test this for both r. and c., since r. descriptions on a given transcript
        have different internal coordinates the  c.
        """
        SDHD = ["ENST00000375549.8", "NC_000011.10(NM_003002.4)", "NM_003002.4"]

        # First, we compare ENST to NC(NM)
        enst, nc_nm = SDHD[:2]
        hgvs_enst = f"{enst}:{coordinate}.{variant}"
        hgvs_nc_nm = f"{nc_nm}:{coordinate}.{variant}"

        cmp_enst = compare_to_wildtype(hgvs_enst)
        cmp_nm = compare_to_wildtype(hgvs_nc_nm)

        assert cmp_enst == cmp_nm

        # Next, we compare ENST to the bare NM
        enst, nm = SDHD[0], SDHD[-1]
        hgvs_enst = f"{enst}:{coordinate}.{variant}"
        hgvs_nm = f"{nm}:{coordinate}.{variant}"

        cmp_enst = compare_to_wildtype(hgvs_enst)
        cmp_nm = compare_to_wildtype(hgvs_nm)

        assert cmp_enst == cmp_nm

    @pytest.mark.parametrize(*SDHD_transcript)
    @pytest.mark.parametrize("variant", SDHD_variants)
    def test_different_coordinates_SDHD(self, transcript: str, variant: str) -> None:
        """Test that the c. and r. variants are handled the same, for ENST, NC(NM) and NM

        Because the internal coordinate systems differ we cannot compare these
        directly but we have to compare them via the comparison with the wildtype
        version of the corresponding transcript

        For each transcript, the r. and c. should give the same predicted
        effect relative to the wildtype.
        """

        hgvs_c = f"{transcript}:c.{variant}"
        hgvs_r = f"{transcript}:r.{variant}"

        cmp_c = compare_to_wildtype(hgvs_c)
        cmp_r = compare_to_wildtype(hgvs_r)

        assert cmp_r == cmp_c

    @pytest.mark.slow
    @pytest.mark.parametrize(*coordinate)
    @pytest.mark.parametrize("variant", WT1_variants)
    def test_different_transcripts_WT1(self, coordinate: str, variant: str) -> None:
        """Test that ENST, NC(NM) and NM variants are handled the same

        Because the internal coordinate systems differ we cannot compare these
        directly but we have to compare them via the comparison with the wildtype
        version of the corresponding transcript

        Test this for both r. and c., since r. descriptions on a given transcript
        have different internal coordinates the  c.
        """
        WT1 = ["ENST00000452863.10", "NC_000011.test(NM_024426.6)", "NM_024426.6"]

        # First, we compare ENST to NC(NM)
        enst, nc_nm = WT1[:2]
        hgvs_enst = f"{enst}:{coordinate}.{variant}"
        hgvs_nc_nm = f"{nc_nm}:{coordinate}.{variant}"

        cmp_enst = compare_to_wildtype(hgvs_enst)
        cmp_nc_nm = compare_to_wildtype(hgvs_nc_nm)

        assert cmp_enst == cmp_nc_nm

        # Next, we compare ENST to the bare NM
        enst, nm = WT1[0], WT1[-1]
        hgvs_enst = f"{enst}:{coordinate}.{variant}"
        hgvs_nm = f"{nm}:{coordinate}.{variant}"

        cmp_enst = compare_to_wildtype(hgvs_enst)
        cmp_nm = compare_to_wildtype(hgvs_nm)

        assert cmp_enst == cmp_nm

    @pytest.mark.parametrize(*WT1_transcript)
    @pytest.mark.parametrize("variant", WT1_variants)
    def test_different_coordinates_WT1(self, transcript: str, variant: str) -> None:
        """Test that the c. and r. variants are handled the same, for ENST, NC(NM) and NM

        Because the internal coordinate systems differ we cannot compare these
        directly but we have to compare them via the comparison with the wildtype
        version of the corresponding transcript

        For each transcript, the r. and c. should give the same predicted
        effect relative to the wildtype.
        """

        hgvs_c = f"{transcript}:c.{variant}"
        hgvs_r = f"{transcript}:r.{variant}"

        cmp_c = compare_to_wildtype(hgvs_c)
        cmp_r = compare_to_wildtype(hgvs_r)

        assert cmp_r == cmp_c
