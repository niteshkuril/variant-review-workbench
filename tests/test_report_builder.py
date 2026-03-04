from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from src.models import (
    AnnotatedVariant,
    ClinVarMatch,
    ConflictSummary,
    DataProvenance,
    GenomeAssembly,
    InputVariant,
    RankedVariant,
    ReviewPriorityTier,
    RunMetadata,
)
from src.report_builder import build_report_context, build_variant_export_records, render_html_report, write_html_report


def build_ranked_variant(
    record_id: str,
    gene: str,
    position: int,
    tier: ReviewPriorityTier,
    score: float,
    *,
    conflict: bool = False,
    matched: bool = True,
) -> RankedVariant:
    input_variant = InputVariant(
        record_id=record_id,
        assembly=GenomeAssembly.GRCH38,
        chromosome="17",
        position=position,
        reference_allele="A",
        alternate_allele="G",
        gene=gene,
        impact="HIGH",
    )
    clinvar = ClinVarMatch(
        matched=matched,
        assembly=GenomeAssembly.GRCH38 if matched else GenomeAssembly.UNKNOWN,
        chromosome="17" if matched else None,
        position=position if matched else None,
        reference_allele="A" if matched else None,
        alternate_allele="G" if matched else None,
        variation_id=position if matched else None,
        gene=gene if matched else None,
        clinical_significance="Pathogenic" if matched else None,
        review_status="reviewed by expert panel" if matched else None,
        review_stars=3 if matched else None,
        condition_names=["Li-Fraumeni syndrome"] if matched else [],
        conflict=ConflictSummary(has_conflict=conflict),
    )
    annotated = AnnotatedVariant(
        input_variant=input_variant,
        clinvar=clinvar,
        flags=["clinvar_matched"] if matched else ["clinvar_unmatched"],
    )
    return RankedVariant(
        annotated_variant=annotated,
        priority_score=score,
        priority_tier=tier,
        ranking_rationale=["Example rationale for report testing."],
    )


class ReportBuilderTests(unittest.TestCase):
    def test_build_variant_export_records_shapes_machine_readable_fields(self) -> None:
        ranked_variants = [
            build_ranked_variant("record-1", "TP53", 43045702, ReviewPriorityTier.HIGH_REVIEW_PRIORITY, 17.5, conflict=True),
        ]

        records = build_variant_export_records(ranked_variants)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].record_id, "record-1")
        self.assertEqual(records[0].input_gene, "TP53")
        self.assertTrue(records[0].clinvar_matched)
        self.assertTrue(records[0].conflict_flagged)
        self.assertEqual(records[0].priority_tier, ReviewPriorityTier.HIGH_REVIEW_PRIORITY)

    def test_build_report_context_summarizes_ranked_variants(self) -> None:
        ranked_variants = [
            build_ranked_variant("record-1", "TP53", 43045702, ReviewPriorityTier.HIGH_REVIEW_PRIORITY, 17.5, conflict=True),
            build_ranked_variant("record-2", "BRAF", 140453136, ReviewPriorityTier.REVIEW, 7.0),
            build_ranked_variant("record-3", "NA", 555, ReviewPriorityTier.CONTEXT_ONLY, 0.0, matched=False),
        ]
        metadata = RunMetadata(
            input_path="data/demo.vcf",
            output_dir="outputs/demo",
            assembly=GenomeAssembly.GRCH38,
            sources=[DataProvenance(source_name="ClinVar variant summary", source_kind="file", source_path="data/clinvar/raw/variant_summary.txt.gz")],
        )

        context = build_report_context(ranked_variants, metadata)

        self.assertEqual(context["summary"]["variant_count"], 3)
        self.assertEqual(context["summary"]["clinvar_matched_count"], 2)
        self.assertEqual(context["summary"]["conflict_count"], 1)
        self.assertEqual(context["summary_artifact"]["input_variant_count"], 3)
        self.assertEqual(context["summary_artifact"]["conflict_flagged_count"], 1)
        self.assertEqual(len(context["top_findings"]), 3)
        self.assertEqual(len(context["conflict_rows"]), 1)
        self.assertEqual(context["variant_rows"][0]["gene"], "TP53")

    def test_render_html_report_contains_expected_sections(self) -> None:
        ranked_variants = [
            build_ranked_variant("record-1", "TP53", 43045702, ReviewPriorityTier.HIGH_REVIEW_PRIORITY, 17.5, conflict=True),
            build_ranked_variant("record-2", "BRAF", 140453136, ReviewPriorityTier.REVIEW, 7.0),
        ]
        metadata = RunMetadata(
            input_path="data/demo.vcf",
            output_dir="outputs/demo",
            assembly=GenomeAssembly.GRCH38,
        )

        html = render_html_report(ranked_variants, metadata)

        self.assertIn("Variant Review Report", html)
        self.assertIn("Top Findings", html)
        self.assertIn("Conflict Review Queue", html)
        self.assertIn("Variant Table", html)
        self.assertIn("Methods", html)
        self.assertIn("Limitations", html)
        self.assertIn("TP53", html)
        self.assertIn("Li-Fraumeni syndrome", html)

    def test_write_html_report_writes_output_file(self) -> None:
        ranked_variants = [
            build_ranked_variant("record-1", "TP53", 43045702, ReviewPriorityTier.HIGH_REVIEW_PRIORITY, 17.5),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "report.html"
            written = write_html_report(output_path, ranked_variants)

            self.assertEqual(written, output_path)
            self.assertTrue(output_path.exists())
            html = output_path.read_text(encoding="utf-8")

        self.assertIn("Variant Review Report", html)
        self.assertIn("TP53", html)


if __name__ == "__main__":
    unittest.main()
