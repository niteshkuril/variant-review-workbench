from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from src.clinvar_index import (
    ClinVarIndex,
    enrich_index_with_supporting_data,
    load_conflict_lookup,
    load_submission_lookup,
    load_variant_summary_index,
    review_status_to_stars,
)
from src.models import ClinVarMatch, ConflictSummary, GenomeAssembly, InputVariant
from src.vcf_parser import normalize_allele, normalize_chromosome, parse_info_field


class ModelAndHelperTests(unittest.TestCase):
    def test_input_variant_key_uses_normalized_shape(self) -> None:
        variant = InputVariant(
            record_id="record-1",
            assembly=GenomeAssembly.GRCH38,
            chromosome="17",
            position=43045702,
            reference_allele="A",
            alternate_allele="G",
        )

        self.assertEqual(variant.variant_key, ("GRCh38", "17", 43045702, "A", "G"))

    def test_clinvar_match_key_defaults_to_unknown_when_unmatched(self) -> None:
        match = ClinVarMatch()

        self.assertEqual(match.variant_key, ("unknown", None, None, None, None))

    def test_parse_info_field_preserves_flags_and_key_values(self) -> None:
        parsed = parse_info_field("GENE=TP53;DB;IMPACT=HIGH")

        self.assertEqual(parsed["GENE"], "TP53")
        self.assertEqual(parsed["DB"], "true")
        self.assertEqual(parsed["IMPACT"], "HIGH")

    def test_normalize_chromosome_and_allele_helpers(self) -> None:
        self.assertEqual(normalize_chromosome("chrM"), "MT")
        self.assertEqual(normalize_chromosome("chr7"), "7")
        self.assertEqual(normalize_allele(" at "), "AT")

    def test_review_status_to_stars_maps_known_values(self) -> None:
        self.assertEqual(review_status_to_stars("practice guideline"), 4)
        self.assertEqual(review_status_to_stars("reviewed by expert panel"), 3)
        self.assertEqual(review_status_to_stars("criteria provided, multiple submitters, no conflicts"), 2)
        self.assertEqual(review_status_to_stars("criteria provided, single submitter"), 1)
        self.assertEqual(review_status_to_stars("no assertion criteria provided"), 0)
        self.assertIsNone(review_status_to_stars(None))

    def test_load_conflict_lookup_filters_to_target_variation_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "conflicts.txt"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(
                    [
                        "#Gene_Symbol",
                        "NCBI_Variation_ID",
                        "ClinVar_Preferred",
                        "Submitter1",
                        "Submitter1_SCV",
                        "Submitter1_ClinSig",
                        "Submitter1_LastEval",
                        "Submitter1_ReviewStatus",
                        "Submitter1_Sub_Condition",
                        "Submitter1_Description",
                        "Submitter2",
                        "Submitter2_SCV",
                        "Submitter2_ClinSig",
                        "Submitter2_LastEval",
                        "Submitter2_ReviewStatus",
                        "Submitter2_Sub_Condition",
                        "Submitter2_Description",
                        "Rank_diff",
                        "Conflict_Reported",
                        "Variant_type",
                        "Submitter1_Method",
                        "Submitter2_Method",
                    ]
                )
                writer.writerow(["TP53", "1", "var1", "A", "SCV1", "Pathogenic", "", "", "", "", "B", "SCV2", "Benign", "", "", "", "", "1", "yes", "SNV", "", ""])
                writer.writerow(["TP53", "2", "var2", "A", "SCV1", "Pathogenic", "", "", "", "", "B", "SCV2", "Benign", "", "", "", "", "1", "yes", "SNV", "", ""])

            conflicts, _ = load_conflict_lookup(path, {2})

        self.assertEqual(set(conflicts.keys()), {2})
        self.assertTrue(conflicts[2].has_conflict)

    def test_load_submission_lookup_skips_preamble_and_filters_targets(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "submission_summary.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                handle.write("##Overview\n")
                handle.write("#VariationID\tClinicalSignificance\tDateLastEvaluated\tDescription\tSubmittedPhenotypeInfo\tReportedPhenotypeInfo\tReviewStatus\tCollectionMethod\tOriginCounts\tSubmitter\tSCV\tSubmittedGeneSymbol\tExplanationOfInterpretation\tSomaticClinicalImpact\tOncogenicity\tContributesToAggregateClassification\n")
                handle.write("1\tPathogenic\tJan 01, 2025\t-\t-\t-\treviewed by expert panel\tclinical testing\tgermline:1\tLab A\tSCV1\tTP53\t-\t-\t-\tyes\n")
                handle.write("2\tBenign\tJan 01, 2025\t-\t-\t-\tcriteria provided, single submitter\tclinical testing\tgermline:1\tLab B\tSCV2\tTP53\t-\t-\t-\tyes\n")

            submissions, _ = load_submission_lookup(path, {1})

        self.assertEqual(set(submissions.keys()), {1})
        self.assertEqual(submissions[1].submitter_names, ["Lab A"])

    def test_enrich_index_with_supporting_data_preserves_base_match_and_copies_on_lookup(self) -> None:
        base_match = ClinVarMatch(
            matched=True,
            assembly=GenomeAssembly.GRCH38,
            chromosome="17",
            position=43045702,
            reference_allele="A",
            alternate_allele="G",
            variation_id=1234,
            review_stars=2,
        )
        index = ClinVarIndex(
            exact_matches={("GRCh38", "17", 43045702, "A", "G"): base_match},
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            conflict_path = Path(tmpdir) / "conflicts.txt"
            with conflict_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(
                    [
                        "#Gene_Symbol",
                        "NCBI_Variation_ID",
                        "ClinVar_Preferred",
                        "Submitter1",
                        "Submitter1_SCV",
                        "Submitter1_ClinSig",
                        "Submitter1_LastEval",
                        "Submitter1_ReviewStatus",
                        "Submitter1_Sub_Condition",
                        "Submitter1_Description",
                        "Submitter2",
                        "Submitter2_SCV",
                        "Submitter2_ClinSig",
                        "Submitter2_LastEval",
                        "Submitter2_ReviewStatus",
                        "Submitter2_Sub_Condition",
                        "Submitter2_Description",
                        "Rank_diff",
                        "Conflict_Reported",
                        "Variant_type",
                        "Submitter1_Method",
                        "Submitter2_Method",
                    ]
                )
                writer.writerow(["TP53", "1234", "var", "A", "SCV1", "Pathogenic", "", "", "", "", "B", "SCV2", "Benign", "", "", "", "", "1", "yes", "SNV", "", ""])

            enrich_index_with_supporting_data(index, conflict_summary_path=conflict_path)

        looked_up = index.lookup(
            InputVariant(
                record_id="record-1",
                assembly=GenomeAssembly.GRCH38,
                chromosome="17",
                position=43045702,
                reference_allele="A",
                alternate_allele="G",
            )
        )
        looked_up.conflict = ConflictSummary(has_conflict=False)

        self.assertTrue(index.conflicts_by_variation_id[1234].has_conflict)

    def test_load_variant_summary_index_skips_malformed_numeric_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "variant_summary.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(
                    [
                        "#AlleleID",
                        "Type",
                        "Name",
                        "GeneSymbol",
                        "ClinicalSignificance",
                        "LastEvaluated",
                        "PhenotypeList",
                        "ReviewStatus",
                        "Origin",
                        "Assembly",
                        "Chromosome",
                        "VariationID",
                        "PositionVCF",
                        "ReferenceAlleleVCF",
                        "AlternateAlleleVCF",
                        "RCVaccession",
                    ]
                )
                writer.writerow(
                    [
                        "",
                        "single nucleotide variant",
                        "bad row",
                        "TP53",
                        "Pathogenic",
                        "Jan 01, 2025",
                        "-",
                        "reviewed by expert panel",
                        "germline",
                        "GRCh38",
                        "17",
                        "",
                        "43045702",
                        "A",
                        "G",
                        "RCV000000001",
                    ]
                )
                writer.writerow(
                    [
                        "10",
                        "single nucleotide variant",
                        "good row",
                        "TP53",
                        "Pathogenic",
                        "Jan 01, 2025",
                        "Li-Fraumeni syndrome",
                        "reviewed by expert panel",
                        "germline",
                        "GRCh38",
                        "17",
                        "1234",
                        "43045702",
                        "A",
                        "G",
                        "RCV000000002",
                    ]
                )

            index = load_variant_summary_index(path)

        self.assertEqual(len(index.exact_matches), 1)
        self.assertIn(("GRCh38", "17", 43045702, "A", "G"), index.exact_matches)


if __name__ == "__main__":
    unittest.main()
