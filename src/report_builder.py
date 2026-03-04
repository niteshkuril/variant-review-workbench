"""HTML and tabular report generation."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from .models import (
    RankedVariant,
    ReviewPriorityTier,
    RunMetadata,
    SummaryArtifact,
    VariantExportRecord,
)

TEMPLATE_DIR = Path(__file__).resolve().parent.parent / "templates"
REPORT_TEMPLATE_NAME = "report.html.j2"

METHODS_NOTES = [
    "Input variants are normalized from VCF into one record per alternate allele.",
    "ClinVar matching uses exact assembly-aware coordinate and allele lookup against the local snapshot.",
    "Conflict and submission summaries are attached by ClinVar VariationID when available.",
    "Ranking is heuristic and transparent; it prioritizes significance, review strength, conflicts, and optional PharmGKB context.",
]

LIMITATIONS_NOTES = [
    "This report is for research triage and educational review only, not clinical decision-making.",
    "Exact matching can miss representation differences that require normalization beyond the current v1 pipeline.",
    "Unmatched variants are not interpreted automatically and should be reviewed manually.",
    "ClinVar significance labels and ranking output should be interpreted in the context of the source snapshot date.",
]


def _build_environment() -> Environment:
    """Create the Jinja environment used for HTML report rendering."""
    return Environment(
        loader=FileSystemLoader(str(TEMPLATE_DIR)),
        autoescape=select_autoescape(["html", "xml"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )


def build_report_summary(ranked_variants: list[RankedVariant]) -> SummaryArtifact:
    """Compute analyst-facing summary metrics for a ranked variant set."""
    tier_counts = Counter(ranked.priority_tier.value for ranked in ranked_variants)
    matched_count = sum(1 for ranked in ranked_variants if ranked.annotated_variant.has_clinvar_match)
    conflict_count = sum(1 for ranked in ranked_variants if ranked.annotated_variant.has_conflict)
    pharmgkb_count = sum(
        1
        for ranked in ranked_variants
        if ranked.annotated_variant.pharmgkb is not None and ranked.annotated_variant.pharmgkb.matched
    )

    return SummaryArtifact(
        input_variant_count=len(ranked_variants),
        clinvar_matched_count=matched_count,
        clinvar_unmatched_count=len(ranked_variants) - matched_count,
        conflict_flagged_count=conflict_count,
        pharmgkb_enriched_count=pharmgkb_count,
        review_priority_tier_counts={
            ReviewPriorityTier.HIGH_REVIEW_PRIORITY.value: tier_counts.get(ReviewPriorityTier.HIGH_REVIEW_PRIORITY.value, 0),
            ReviewPriorityTier.REVIEW.value: tier_counts.get(ReviewPriorityTier.REVIEW.value, 0),
            ReviewPriorityTier.CONTEXT_ONLY.value: tier_counts.get(ReviewPriorityTier.CONTEXT_ONLY.value, 0),
        },
    )


def _format_conditions(conditions: list[str]) -> str:
    """Render ClinVar conditions as analyst-friendly text."""
    return ", ".join(conditions) if conditions else "Unspecified"


def _build_variant_row(ranked_variant: RankedVariant) -> dict[str, object]:
    """Convert a RankedVariant into a report/table-friendly record."""
    annotated = ranked_variant.annotated_variant
    input_variant = annotated.input_variant
    clinvar = annotated.clinvar
    pharmgkb = annotated.pharmgkb

    return {
        "record_id": input_variant.record_id,
        "locus": f"{input_variant.chromosome}:{input_variant.position} {input_variant.reference_allele}>{input_variant.alternate_allele}",
        "gene": clinvar.gene or input_variant.gene or "Unspecified",
        "priority_score": ranked_variant.priority_score,
        "priority_tier": ranked_variant.priority_tier.value,
        "clinical_significance": clinvar.clinical_significance or "No ClinVar match",
        "review_status": clinvar.review_status or "Unspecified",
        "review_stars": clinvar.review_stars,
        "condition_names": _format_conditions(clinvar.condition_names),
        "conflict": "Yes" if annotated.has_conflict else "No",
        "variation_id": clinvar.variation_id,
        "preferred_name": clinvar.preferred_name or "Unspecified",
        "impact": input_variant.impact or "Unspecified",
        "transcript": input_variant.transcript or "Unspecified",
        "pharmgkb": "Yes" if pharmgkb is not None and pharmgkb.matched else "No",
        "flags": annotated.flags,
        "rationale": ranked_variant.ranking_rationale,
    }


def build_variant_export_records(ranked_variants: list[RankedVariant]) -> list[VariantExportRecord]:
    """Build stable machine-readable export records from ranked variants."""
    records: list[VariantExportRecord] = []
    for ranked_variant in ranked_variants:
        annotated = ranked_variant.annotated_variant
        input_variant = annotated.input_variant
        clinvar = annotated.clinvar
        pharmgkb = annotated.pharmgkb
        submissions = clinvar.submissions

        records.append(
            VariantExportRecord(
                record_id=input_variant.record_id,
                assembly=input_variant.assembly,
                chromosome=input_variant.chromosome,
                position=input_variant.position,
                reference_allele=input_variant.reference_allele,
                alternate_allele=input_variant.alternate_allele,
                locus=f"{input_variant.chromosome}:{input_variant.position} {input_variant.reference_allele}>{input_variant.alternate_allele}",
                variant_id=input_variant.variant_id,
                input_gene=input_variant.gene,
                clinvar_gene=clinvar.gene,
                preferred_name=clinvar.preferred_name,
                variation_id=clinvar.variation_id,
                allele_id=clinvar.allele_id,
                match_strategy=clinvar.match_strategy,
                clinvar_matched=clinvar.matched,
                clinical_significance=clinvar.clinical_significance,
                review_status=clinvar.review_status,
                review_stars=clinvar.review_stars,
                condition_names=list(clinvar.condition_names),
                conflict_flagged=annotated.has_conflict,
                conflict_significance=list(clinvar.conflict.conflict_significance),
                conflict_summary_text=clinvar.conflict.summary_text,
                submission_count=submissions.total_submissions if submissions is not None else None,
                submission_submitter_names=list(submissions.submitter_names) if submissions is not None else [],
                submission_review_statuses=list(submissions.review_statuses) if submissions is not None else [],
                submission_clinical_significances=list(submissions.clinical_significances) if submissions is not None else [],
                pharmgkb_queried=pharmgkb.queried if pharmgkb is not None else False,
                pharmgkb_matched=pharmgkb.matched if pharmgkb is not None else False,
                pharmgkb_gene_ids=list(pharmgkb.pharmgkb_gene_ids) if pharmgkb is not None else [],
                pharmgkb_variant_ids=list(pharmgkb.pharmgkb_variant_ids) if pharmgkb is not None else [],
                pharmgkb_chemicals=list(pharmgkb.chemicals) if pharmgkb is not None else [],
                input_transcript=input_variant.transcript,
                input_impact=input_variant.impact,
                input_consequence=input_variant.consequence,
                priority_score=ranked_variant.priority_score,
                priority_tier=ranked_variant.priority_tier,
                flags=list(annotated.flags),
                ranking_rationale=list(ranked_variant.ranking_rationale),
            )
        )
    return records


def build_report_context(
    ranked_variants: list[RankedVariant],
    run_metadata: RunMetadata | None = None,
) -> dict[str, object]:
    """Build the full template context for the HTML report."""
    summary = build_report_summary(ranked_variants)
    rows = [_build_variant_row(ranked_variant) for ranked_variant in ranked_variants]
    top_findings = rows[: min(5, len(rows))]
    conflict_rows = [row for row in rows if row["conflict"] == "Yes"]
    report_summary = {
        "variant_count": summary.input_variant_count,
        "clinvar_matched_count": summary.clinvar_matched_count,
        "clinvar_unmatched_count": summary.clinvar_unmatched_count,
        "conflict_count": summary.conflict_flagged_count,
        "pharmgkb_count": summary.pharmgkb_enriched_count,
        "tier_counts": summary.review_priority_tier_counts,
    }

    return {
        "report_title": "Variant Review Report",
        "generated_at": run_metadata.run_started_at.isoformat() if run_metadata is not None else None,
        "assembly": run_metadata.assembly.value if run_metadata is not None else None,
        "input_path": run_metadata.input_path if run_metadata is not None else None,
        "summary": report_summary,
        "summary_artifact": summary.model_dump(mode="json"),
        "top_findings": top_findings,
        "conflict_rows": conflict_rows,
        "variant_rows": rows,
        "methods_notes": METHODS_NOTES,
        "limitations_notes": LIMITATIONS_NOTES,
        "sources": run_metadata.sources if run_metadata is not None else [],
    }


def render_html_report(
    ranked_variants: list[RankedVariant],
    run_metadata: RunMetadata | None = None,
) -> str:
    """Render the analyst-facing HTML report."""
    environment = _build_environment()
    template = environment.get_template(REPORT_TEMPLATE_NAME)
    context = build_report_context(ranked_variants, run_metadata=run_metadata)
    return template.render(**context)


def write_html_report(
    output_path: Path,
    ranked_variants: list[RankedVariant],
    run_metadata: RunMetadata | None = None,
) -> Path:
    """Render and write an HTML report to disk."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        render_html_report(ranked_variants, run_metadata=run_metadata),
        encoding="utf-8",
    )
    return output_path
