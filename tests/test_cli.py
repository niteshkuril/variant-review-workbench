from __future__ import annotations

import argparse
import csv
import gzip
import json
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from src.cli import build_parser, main, run_pipeline
from src.models import RunMetadata
from src.models import GenomeAssembly


class CliTests(unittest.TestCase):
    def test_build_parser_accepts_supported_assembly(self) -> None:
        parser = build_parser()

        args = parser.parse_args(
            [
                "--input",
                "input.vcf",
                "--assembly",
                "GRCh38",
                "--variant-summary",
                "variant_summary.txt.gz",
                "--out-dir",
                "outputs",
            ]
        )

        self.assertEqual(args.assembly, GenomeAssembly.GRCH38)

    def test_build_parser_rejects_unsupported_assembly(self) -> None:
        parser = build_parser()

        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--input",
                    "input.vcf",
                    "--assembly",
                    "hg19",
                    "--variant-summary",
                    "variant_summary.txt.gz",
                    "--out-dir",
                    "outputs",
                ]
            )

    def test_run_pipeline_writes_expected_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input.vcf"
            input_path.write_text(
                (
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "17\t43045702\t.\tA\tG\t100\tPASS\tGENE=TP53;IMPACT=HIGH\n"
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

            out_dir = root / "outputs"
            args = argparse.Namespace(
                input=str(input_path),
                assembly=GenomeAssembly.GRCH38,
                variant_summary=str(variant_summary),
                conflict_summary=str(conflict),
                submission_summary=str(submission),
                clinvar_cache_db=None,
                disable_clinvar_cache=False,
                out_dir=str(out_dir),
                enable_pharmgkb=False,
            )

            outputs = run_pipeline(args)

            for path in outputs.values():
                self.assertTrue(path.exists())

            summary = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
            metadata = json.loads((out_dir / "run_metadata.json").read_text(encoding="utf-8"))
            report_html = (out_dir / "report.html").read_text(encoding="utf-8")
            variants_json = json.loads((out_dir / "prioritized_variants.json").read_text(encoding="utf-8"))
            with (out_dir / "annotated_variants.csv").open("r", encoding="utf-8", newline="") as handle:
                csv_rows = list(csv.DictReader(handle))

        self.assertEqual(summary["schema_version"], "1.0")
        self.assertEqual(summary["artifact_type"], "summary")
        self.assertEqual(summary["input_variant_count"], 1)
        self.assertEqual(summary["conflict_flagged_count"], 1)
        self.assertEqual(metadata["statistics"]["clinvar_matched_count"], 1)
        self.assertIn("Variant Review Report", report_html)
        self.assertEqual(variants_json["schema_version"], "1.0")
        self.assertEqual(variants_json["artifact_type"], "prioritized_variants")
        self.assertEqual(variants_json["records"][0]["input_gene"], "TP53")
        self.assertTrue(variants_json["records"][0]["conflict_flagged"])
        self.assertEqual(variants_json["records"][0]["condition_names"], ["Li-Fraumeni syndrome"])
        self.assertEqual(csv_rows[0]["clinvar_matched"], "true")
        self.assertEqual(csv_rows[0]["condition_names"], "[\"Li-Fraumeni syndrome\"]")
        self.assertEqual(csv_rows[0]["flags"], "[\"clinvar_matched\", \"clinvar_conflict\", \"clinvar_review_stars_3\", \"submission_evidence_available\"]")

    def test_run_pipeline_raises_clear_error_for_missing_input_file(self) -> None:
        args = argparse.Namespace(
            input="missing.vcf",
            assembly=GenomeAssembly.GRCH38,
            variant_summary="variant_summary.txt.gz",
            conflict_summary=None,
            submission_summary=None,
            clinvar_cache_db=None,
            disable_clinvar_cache=False,
            out_dir="outputs",
            enable_pharmgkb=False,
        )

        with self.assertRaisesRegex(ValueError, "Input VCF was not found: missing.vcf"):
            run_pipeline(args)

    def test_main_reports_when_pharmgkb_enabled_but_no_matches_found(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input.vcf"
            input_path.write_text(
                (
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "16\t31105878\trs104894541\tC\tT\t100\tPASS\tGENE=VKORC1;IMPACT=HIGH\n"
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
                        "VKORC1 example",
                        "VKORC1",
                        "drug response",
                        "Jan 01, 2025",
                        "Warfarin response",
                        "criteria provided, single submitter",
                        "germline",
                        "GRCh38",
                        "16",
                        "1234",
                        "31105878",
                        "C",
                        "T",
                        "RCV000000001",
                    ]
                )

            out_dir = root / "outputs"
            argv = [
                "prog",
                "--input",
                str(input_path),
                "--assembly",
                "GRCh38",
                "--variant-summary",
                str(variant_summary),
                "--out-dir",
                str(out_dir),
                "--enable-pharmgkb",
            ]

            with patch("src.pgx_enrichment.PharmGKBClient._get", return_value=([], False)):
                with patch("sys.argv", argv):
                    with patch("sys.stdout", new_callable=io.StringIO) as stdout:
                        main()

        output = stdout.getvalue()
        self.assertIn("Run completed:", output)
        self.assertIn("PharmGKB was enabled but no enrichment matches were found.", output)

    def test_cli_main_uses_shared_pipeline_service(self) -> None:
        outputs = {"report_html": Path("report.html")}
        run_metadata = RunMetadata(
            input_path="input.vcf",
            output_dir="outputs",
            assembly=GenomeAssembly.GRCH38,
            pharmgkb_enabled=False,
            sources=[],
        )
        run_metadata.statistics.input_variant_count = 1
        run_metadata.statistics.clinvar_matched_count = 1
        run_metadata.statistics.conflict_flagged_count = 0
        run_metadata.statistics.pharmgkb_enriched_count = 0

        with patch("src.cli.run_pipeline_with_details", return_value=(outputs, run_metadata)) as mock_run:
            with patch("sys.argv", ["prog", "--input", "input.vcf", "--assembly", "GRCh38", "--variant-summary", "variant_summary.txt.gz", "--out-dir", "outputs"]):
                with patch("sys.stdout", new_callable=io.StringIO) as stdout:
                    main()

        self.assertTrue(mock_run.called)
        self.assertIn("Run completed:", stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
