from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from src.annotator import annotate_variant, annotate_variants
from src.clinvar_index import load_clinvar_index
from src.models import GenomeAssembly
from src.vcf_parser import parse_vcf


class AnnotatorTests(unittest.TestCase):
    def test_annotate_variant_marks_unmatched_variants(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input.vcf"
            input_path.write_text(
                (
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "17\t43045702\t.\tA\tG\t100\tPASS\tGENE=TP53\n"
                ),
                encoding="utf-8",
            )

            variant_summary = root / "variant_summary.txt.gz"
            with gzip.open(variant_summary, "wt", encoding="utf-8", newline="") as handle:
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

            index = load_clinvar_index(variant_summary)
            variant = parse_vcf(input_path, GenomeAssembly.GRCH38)[0]
            annotated = annotate_variant(variant, index)

        self.assertFalse(annotated.has_clinvar_match)
        self.assertEqual(annotated.flags, ["clinvar_unmatched"])

    def test_annotate_variant_attaches_conflict_and_submission_flags(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input.vcf"
            input_path.write_text(
                (
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "17\t43045702\t.\tA\tG\t100\tPASS\tGENE=TP53\n"
                ),
                encoding="utf-8",
            )

            variant_summary = root / "variant_summary.txt.gz"
            with gzip.open(variant_summary, "wt", encoding="utf-8", newline="") as handle:
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
                        "10",
                        "single nucleotide variant",
                        "TP53 example",
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
                        "RCV000000001",
                    ]
                )

            conflict = root / "summary_of_conflicting_interpretations.txt"
            with conflict.open("w", encoding="utf-8", newline="") as handle:
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
                writer.writerow(
                    [
                        "TP53",
                        "1234",
                        "TP53 example",
                        "Lab A",
                        "SCV1",
                        "Pathogenic",
                        "Jan 01, 2025",
                        "criteria provided, single submitter",
                        "",
                        "",
                        "Lab B",
                        "SCV2",
                        "Uncertain significance",
                        "Jan 02, 2025",
                        "criteria provided, single submitter",
                        "",
                        "",
                        "1",
                        "yes",
                        "SNV",
                        "clinical testing",
                        "clinical testing",
                    ]
                )

            submission = root / "submission_summary.txt.gz"
            with gzip.open(submission, "wt", encoding="utf-8", newline="") as handle:
                handle.write("##Comment\n")
                handle.write(
                    "#VariationID\tClinicalSignificance\tDateLastEvaluated\tDescription\t"
                    "SubmittedPhenotypeInfo\tReportedPhenotypeInfo\tReviewStatus\tCollectionMethod\t"
                    "OriginCounts\tSubmitter\tSCV\tSubmittedGeneSymbol\tExplanationOfInterpretation\t"
                    "SomaticClinicalImpact\tOncogenicity\tContributesToAggregateClassification\n"
                )
                handle.write(
                    "1234\tPathogenic\tJan 01, 2025\t-\tLi-Fraumeni syndrome\t"
                    "C0085390:Li-Fraumeni syndrome\treviewed by expert panel\tclinical testing\t"
                    "germline:1\tLab A\tSCV1\tTP53\t-\t-\t-\tyes\n"
                )

            index = load_clinvar_index(variant_summary, conflict, submission)
            variant = parse_vcf(input_path, GenomeAssembly.GRCH38)[0]
            annotated = annotate_variant(variant, index)

        self.assertTrue(annotated.has_clinvar_match)
        self.assertTrue(annotated.has_conflict)
        self.assertIn("clinvar_matched", annotated.flags)
        self.assertIn("clinvar_conflict", annotated.flags)
        self.assertIn("clinvar_review_stars_3", annotated.flags)
        self.assertIn("submission_evidence_available", annotated.flags)

    def test_annotate_variant_flags_gene_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input.vcf"
            input_path.write_text(
                (
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "17\t43045702\t.\tA\tG\t100\tPASS\tGENE=WRONG1\n"
                ),
                encoding="utf-8",
            )

            variant_summary = root / "variant_summary.txt.gz"
            with gzip.open(variant_summary, "wt", encoding="utf-8", newline="") as handle:
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
                        "10",
                        "single nucleotide variant",
                        "TP53 example",
                        "TP53",
                        "Pathogenic",
                        "Jan 01, 2025",
                        "Li-Fraumeni syndrome",
                        "criteria provided, single submitter",
                        "germline",
                        "GRCh38",
                        "17",
                        "1234",
                        "43045702",
                        "A",
                        "G",
                        "RCV000000001",
                    ]
                )

            index = load_clinvar_index(variant_summary)
            variant = parse_vcf(input_path, GenomeAssembly.GRCH38)[0]
            annotated = annotate_variant(variant, index)

        self.assertTrue(annotated.has_clinvar_match)
        self.assertIn("gene_symbol_mismatch", annotated.flags)

    def test_annotate_variant_does_not_flag_case_only_gene_difference(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input.vcf"
            input_path.write_text(
                (
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "17\t43045702\t.\tA\tG\t100\tPASS\tGENE=tp53\n"
                ),
                encoding="utf-8",
            )

            variant_summary = root / "variant_summary.txt.gz"
            with gzip.open(variant_summary, "wt", encoding="utf-8", newline="") as handle:
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
                        "10",
                        "single nucleotide variant",
                        "TP53 example",
                        "TP53",
                        "Pathogenic",
                        "Jan 01, 2025",
                        "Li-Fraumeni syndrome",
                        "criteria provided, single submitter",
                        "germline",
                        "GRCh38",
                        "17",
                        "1234",
                        "43045702",
                        "A",
                        "G",
                        "RCV000000001",
                    ]
                )

            index = load_clinvar_index(variant_summary)
            variant = parse_vcf(input_path, GenomeAssembly.GRCH38)[0]
            annotated = annotate_variant(variant, index)

        self.assertTrue(annotated.has_clinvar_match)
        self.assertNotIn("gene_symbol_mismatch", annotated.flags)

    def test_annotate_variants_preserves_input_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input.vcf"
            input_path.write_text(
                (
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "17\t43045702\t.\tA\tG\t100\tPASS\tGENE=TP53\n"
                    "7\t140453136\t.\tA\tT\t100\tPASS\tGENE=BRAF\n"
                ),
                encoding="utf-8",
            )

            variant_summary = root / "variant_summary.txt.gz"
            with gzip.open(variant_summary, "wt", encoding="utf-8", newline="") as handle:
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
                        "10",
                        "single nucleotide variant",
                        "TP53 example",
                        "TP53",
                        "Pathogenic",
                        "Jan 01, 2025",
                        "Li-Fraumeni syndrome",
                        "criteria provided, single submitter",
                        "germline",
                        "GRCh38",
                        "17",
                        "1234",
                        "43045702",
                        "A",
                        "G",
                        "RCV000000001",
                    ]
                )

            index = load_clinvar_index(variant_summary)
            variants = parse_vcf(input_path, GenomeAssembly.GRCH38)
            annotated_variants = annotate_variants(variants, index)

        self.assertEqual(
            [item.input_variant.record_id for item in annotated_variants],
            ["record-1", "record-2"],
        )
        self.assertTrue(annotated_variants[0].has_clinvar_match)
        self.assertFalse(annotated_variants[1].has_clinvar_match)


if __name__ == "__main__":
    unittest.main()
