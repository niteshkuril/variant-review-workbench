"""HTML and tabular report generation."""

from __future__ import annotations

from collections import Counter
import json
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
MAX_REPORT_TABLE_ROWS = 500

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


def _build_no_clinvar_match_warning(summary: SummaryArtifact, assembly: str | None) -> dict[str, str] | None:
    """Build a user-facing warning when no variants matched ClinVar."""
    if summary.input_variant_count == 0 or summary.clinvar_matched_count > 0:
        return None

    assembly_label = assembly or "selected assembly"
    return {
        "title": "No ClinVar matches were found for this run.",
        "message": (
            f"All variants were assigned context_only. Confirm the selected assembly ({assembly_label}) matches the "
            "VCF coordinates (for example, GRCh37 versus GRCh38), and verify the local ClinVar snapshot covers "
            "these loci. Custom assemblies can only be run via command line, requiring a cold start."
        ),
    }


def _serialize_sources(sources: list[object]) -> list[dict[str, object]]:
    """Convert provenance models or mapping-like objects into JSON-safe dicts."""
    serialized: list[dict[str, object]] = []
    for source in sources:
        if hasattr(source, "model_dump"):
            serialized.append(source.model_dump(mode="json"))
        elif isinstance(source, dict):
            serialized.append(dict(source))
    return serialized


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
    assembly_value = run_metadata.assembly.value if run_metadata is not None else None
    rows = [_build_variant_row(ranked_variant) for ranked_variant in ranked_variants[:MAX_REPORT_TABLE_ROWS]]
    top_findings = rows[: min(5, len(rows))]
    conflict_rows = [row for row in rows if row["conflict"] == "Yes"]
    truncated_variant_count = max(0, len(ranked_variants) - MAX_REPORT_TABLE_ROWS)
    report_summary = {
        "variant_count": summary.input_variant_count,
        "displayed_variant_count": len(rows),
        "truncated_variant_count": truncated_variant_count,
        "clinvar_matched_count": summary.clinvar_matched_count,
        "clinvar_unmatched_count": summary.clinvar_unmatched_count,
        "conflict_count": summary.conflict_flagged_count,
        "pharmgkb_count": summary.pharmgkb_enriched_count,
        "tier_counts": summary.review_priority_tier_counts,
    }

    return {
        "report_title": "Variant Review Report",
        "generated_at": run_metadata.run_started_at.isoformat() if run_metadata is not None else None,
        "assembly": assembly_value,
        "input_path": run_metadata.input_path if run_metadata is not None else None,
        "summary": report_summary,
        "summary_artifact": summary.model_dump(mode="json"),
        "no_clinvar_match_warning": _build_no_clinvar_match_warning(summary, assembly_value),
        "top_findings": top_findings,
        "conflict_rows": conflict_rows,
        "variant_rows": rows,
        "methods_notes": METHODS_NOTES,
        "limitations_notes": LIMITATIONS_NOTES,
        "sources": run_metadata.sources if run_metadata is not None else [],
    }


def build_report_export_payload(report_context: dict[str, object]) -> dict[str, object]:
    """Build a JSON-safe report export payload from shared report context."""
    return {
        "report_title": report_context["report_title"],
        "generated_at": report_context["generated_at"],
        "assembly": report_context["assembly"],
        "input_path": report_context["input_path"],
        "summary": report_context["summary"],
        "summary_artifact": report_context["summary_artifact"],
        "no_clinvar_match_warning": report_context.get("no_clinvar_match_warning"),
        "top_findings": report_context["top_findings"],
        "conflict_rows": report_context["conflict_rows"],
        "variant_rows": report_context["variant_rows"],
        "methods_notes": report_context["methods_notes"],
        "limitations_notes": report_context["limitations_notes"],
        "sources": _serialize_sources(list(report_context.get("sources", []))),
    }


def render_markdown_report_from_context(report_context: dict[str, object]) -> str:
    """Render a Markdown report from shared report context."""
    payload = build_report_export_payload(report_context)
    summary = payload["summary"]
    assert isinstance(summary, dict)
    lines = [
        f"# {payload['report_title']}",
        "",
        f"- Generated: {payload['generated_at'] or 'Unspecified'}",
        f"- Assembly: {payload['assembly'] or 'Unspecified'}",
        f"- Input: {payload['input_path'] or 'Unspecified'}",
        "",
        "## Summary",
        "",
        f"- Total variants: {summary['variant_count']}",
        f"- ClinVar matched: {summary['clinvar_matched_count']}",
        f"- ClinVar unmatched: {summary['clinvar_unmatched_count']}",
        f"- Conflict flagged: {summary['conflict_count']}",
        f"- PharmGKB enriched: {summary['pharmgkb_count']}",
        "",
    ]
    warning = payload.get("no_clinvar_match_warning")
    if warning:
        assert isinstance(warning, dict)
        lines.extend(
            [
                "## Match Warning",
                "",
                f"- {warning.get('title', 'No ClinVar matches were found for this run.')}",
                f"- {warning.get('message', '')}",
                "",
            ]
        )
    lines.extend(
        [
        "## Top Findings",
        "",
        ]
    )
    top_findings = payload["top_findings"]
    assert isinstance(top_findings, list)
    for row in top_findings:
        assert isinstance(row, dict)
        lines.extend(
            [
                f"### {row['gene']} ({row['locus']})",
                f"- Priority tier: {row['priority_tier']}",
                f"- Priority score: {row['priority_score']}",
                f"- ClinVar significance: {row['clinical_significance']}",
                f"- Review status: {row['review_status']}",
                f"- Conditions: {row['condition_names']}",
                f"- Conflict: {row['conflict']}",
                "",
            ]
        )

    lines.extend(["## Methods", ""])
    methods_notes = payload["methods_notes"]
    assert isinstance(methods_notes, list)
    lines.extend(f"- {note}" for note in methods_notes)
    lines.extend(["", "## Limitations", ""])
    limitation_notes = payload["limitations_notes"]
    assert isinstance(limitation_notes, list)
    lines.extend(f"- {note}" for note in limitation_notes)
    lines.extend(["", "## Sources", ""])
    sources = payload["sources"]
    assert isinstance(sources, list)
    if sources:
        for source in sources:
            assert isinstance(source, dict)
            location = source.get("source_path") or source.get("source_url") or "Unspecified"
            lines.append(f"- {source.get('source_name', 'Source')}: {location}")
    else:
        lines.append("- No source metadata recorded.")
    lines.append("")
    return "\n".join(lines)


def render_markdown_report(
    ranked_variants: list[RankedVariant],
    run_metadata: RunMetadata | None = None,
) -> str:
    """Render a Markdown report from ranked variants."""
    return render_markdown_report_from_context(build_report_context(ranked_variants, run_metadata=run_metadata))


def write_markdown_report(output_path: Path, report_context: dict[str, object]) -> Path:
    """Write a Markdown report to disk from shared report context."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_markdown_report_from_context(report_context), encoding="utf-8")
    return output_path


def write_report_export_json(output_path: Path, report_context: dict[str, object]) -> Path:
    """Write a JSON report export to disk from shared report context."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(build_report_export_payload(report_context), indent=2), encoding="utf-8")
    return output_path


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
