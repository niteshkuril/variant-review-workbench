from __future__ import annotations

import argparse
import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path

from src.cli import run_pipeline
from src.models import GenomeAssembly


class CliTests(unittest.TestCase):
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

        self.assertEqual(summary["variant_count"], 1)
        self.assertEqual(summary["conflict_count"], 1)
        self.assertEqual(metadata["statistics"]["clinvar_matched_count"], 1)
        self.assertIn("Variant Review Report", report_html)
        self.assertEqual(variants_json[0]["gene"], "TP53")
        self.assertEqual(variants_json[0]["conflict"], "Yes")


if __name__ == "__main__":
    unittest.main()
