from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from src.clinvar_index import load_clinvar_index
from src.models import GenomeAssembly
from src.vcf_parser import parse_vcf


class FoundationPipelineTests(unittest.TestCase):
    def test_parse_vcf_rejects_missing_header(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.vcf"
            input_path.write_text(
                "17\t43045702\t.\tA\tG\t100\tPASS\tGENE=TP53\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "#CHROM header"):
                parse_vcf(input_path, GenomeAssembly.GRCH38)

    def test_parse_vcf_rejects_short_record(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.vcf"
            input_path.write_text(
                (
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "17\t43045702\t.\tA\tG\t100\tPASS\n"
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "fewer than 8 columns"):
                parse_vcf(input_path, GenomeAssembly.GRCH38)

    def test_parse_vcf_normalizes_chromosome_and_allele_case(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.vcf"
            input_path.write_text(
                (
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "chr17\t43045702\t.\ta\tg\t100\tPASS\tGENE=TP53;IMPACT=HIGH\n"
                ),
                encoding="utf-8",
            )

            variants = parse_vcf(input_path, GenomeAssembly.GRCH38)

        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0].chromosome, "17")
        self.assertEqual(variants[0].reference_allele, "A")
        self.assertEqual(variants[0].alternate_allele, "G")

    def test_parse_vcf_expands_multiallelic_records_and_extracts_ann_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.vcf.gz"
            with gzip.open(input_path, "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                handle.write(
                    (
                        "chr7\t140453136\trs123\tc\tt,g\t99\tPASS\t"
                        "ANN=T|missense_variant|MODERATE|BRAF|Gene|transcript|NM_004333.6\n"
                    )
                )

            variants = parse_vcf(input_path, GenomeAssembly.GRCH38)

        self.assertEqual(len(variants), 2)
        self.assertEqual([variant.alternate_allele for variant in variants], ["T", "G"])
        self.assertEqual(variants[0].gene, "BRAF")
        self.assertEqual(variants[0].transcript, "NM_004333.6")
        self.assertEqual(variants[0].consequence, "missense_variant")
        self.assertEqual(variants[0].impact, "MODERATE")
        self.assertEqual(variants[0].variant_id, "rs123")

    def test_parse_vcf_respects_max_variants_limit(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.vcf"
            input_path.write_text(
                (
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "17\t43045702\t.\tA\tG\t100\tPASS\tGENE=TP53;IMPACT=HIGH\n"
                    "17\t43045703\t.\tA\tG\t100\tPASS\tGENE=TP53;IMPACT=HIGH\n"
                    "17\t43045704\t.\tA\tG\t100\tPASS\tGENE=TP53;IMPACT=HIGH\n"
                ),
                encoding="utf-8",
            )

            variants = parse_vcf(input_path, GenomeAssembly.GRCH38, max_variants=2)

        self.assertEqual(len(variants), 2)
        self.assertEqual(variants[0].position, 43045702)
        self.assertEqual(variants[1].position, 43045703)

    def test_clinvar_index_matches_lowercase_input_and_filters_placeholder_conditions(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input.vcf"
            input_path.write_text(
                (
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "chr17\t43045702\t.\ta\tg\t100\tPASS\tGENE=TP53;IMPACT=HIGH\n"
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
                        "not provided|Li-Fraumeni syndrome",
                        "reviewed by expert panel",
                        "germline",
                        "GRCh38",
                        "chr17",
                        "1234",
                        "43045702",
                        "a",
                        "g",
                        "RCV000000001",
                    ]
                )

            index = load_clinvar_index(variant_summary)
            variants = parse_vcf(input_path, GenomeAssembly.GRCH38)
            match = index.lookup(variants[0])

        self.assertTrue(match.matched)
        self.assertEqual(match.reference_allele, "A")
        self.assertEqual(match.alternate_allele, "G")
        self.assertEqual(match.condition_names, ["Li-Fraumeni syndrome"])

    def test_clinvar_index_prefers_stronger_review_status_for_duplicate_keys(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
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
                        "weaker aggregate",
                        "TP53",
                        "Pathogenic",
                        "Jan 01, 2025",
                        "Li-Fraumeni syndrome",
                        "criteria provided, single submitter",
                        "germline",
                        "GRCh38",
                        "17",
                        "111",
                        "43045702",
                        "A",
                        "G",
                        "RCV000000010",
                    ]
                )
                writer.writerow(
                    [
                        "11",
                        "single nucleotide variant",
                        "stronger aggregate",
                        "TP53",
                        "Pathogenic",
                        "Jan 02, 2025",
                        "Li-Fraumeni syndrome",
                        "reviewed by expert panel",
                        "germline",
                        "GRCh38",
                        "17",
                        "222",
                        "43045702",
                        "A",
                        "G",
                        "RCV000000011",
                    ]
                )

            index = load_clinvar_index(variant_summary)
            key = ("GRCh38", "17", 43045702, "A", "G")
            match = index.exact_matches[key]

        self.assertEqual(match.variation_id, 222)
        self.assertEqual(match.review_stars, 3)
        self.assertEqual(match.preferred_name, "stronger aggregate")

    def test_clinvar_index_attaches_conflict_and_submission_support(self) -> None:
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
            variants = parse_vcf(input_path, GenomeAssembly.GRCH38)
            match = index.lookup(variants[0])

        self.assertTrue(match.conflict.has_conflict)
        self.assertEqual(match.conflict.submitter_count, 2)
        self.assertEqual(match.conflict.conflict_significance, ["Pathogenic", "Uncertain significance"])
        self.assertIsNotNone(match.submissions)
        self.assertEqual(match.submissions.total_submissions, 1)
        self.assertEqual(match.submissions.submitter_names, ["Lab A"])


if __name__ == "__main__":
    unittest.main()
