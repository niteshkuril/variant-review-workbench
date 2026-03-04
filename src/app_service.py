"""Shared application service for executing the variant review pipeline."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

from .annotator import annotate_variants
from .clinvar_index import load_clinvar_index
from .models import (
    PrioritizedVariantsArtifact,
    RunMetadata,
    SummaryArtifact,
)
from .pgx_enrichment import PharmGKBClient, enrich_annotated_variants
from .ranker import rank_variants
from .report_builder import build_report_context, build_variant_export_records, write_html_report
from .vcf_parser import parse_vcf


class PipelineUsageError(ValueError):
    """Raised when execution inputs are syntactically valid but operationally unusable."""


@dataclass(slots=True)
class PipelineRunResult:
    """Structured result returned by the shared pipeline service."""

    outputs: dict[str, Path]
    run_metadata: RunMetadata
    report_context: dict[str, object]


def build_run_metadata(args: argparse.Namespace, output_dir: Path, index_sources: list) -> RunMetadata:
    """Construct run metadata after the ClinVar index has been prepared."""
    return RunMetadata(
        input_path=str(Path(args.input)),
        output_dir=str(output_dir),
        assembly=args.assembly,
        pharmgkb_enabled=bool(args.enable_pharmgkb),
        sources=index_sources,
    )


def _validate_existing_file(path: Path, label: str) -> None:
    """Validate that an expected input path exists as a file."""
    if not path.exists():
        raise PipelineUsageError(f"{label} was not found: {path}")
    if not path.is_file():
        raise PipelineUsageError(f"{label} must be a file path: {path}")


def validate_runtime_paths(args: argparse.Namespace) -> None:
    """Validate user-supplied paths before executing the pipeline."""
    _validate_existing_file(Path(args.input), "Input VCF")
    _validate_existing_file(Path(args.variant_summary), "ClinVar variant summary")
    if args.conflict_summary:
        _validate_existing_file(Path(args.conflict_summary), "ClinVar conflict summary")
    if args.submission_summary:
        _validate_existing_file(Path(args.submission_summary), "ClinVar submission summary")

    output_dir = Path(args.out_dir)
    if output_dir.exists() and not output_dir.is_dir():
        raise PipelineUsageError(f"Output path must be a directory: {output_dir}")

    if args.clinvar_cache_db:
        cache_db_path = Path(args.clinvar_cache_db)
        if cache_db_path.exists() and cache_db_path.is_dir():
            raise PipelineUsageError(f"ClinVar cache path must be a file path, not a directory: {cache_db_path}")


def write_json(output_path: Path, payload: object) -> Path:
    """Write a JSON artifact to disk."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return output_path


def write_csv(output_path: Path, rows: list[dict[str, object]]) -> Path:
    """Write a stable CSV export of variant rows."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        output_path.write_text("", encoding="utf-8")
        return output_path

    header = list(rows[0].keys())
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        for row in rows:
            serialized = {}
            for key, value in row.items():
                if isinstance(value, list):
                    serialized[key] = json.dumps(value)
                elif isinstance(value, bool):
                    serialized[key] = "true" if value else "false"
                elif value is None:
                    serialized[key] = ""
                else:
                    serialized[key] = value
            writer.writerow(serialized)
    return output_path


def run_pipeline_with_result(args: argparse.Namespace) -> PipelineRunResult:
    """Execute the local ClinVar-first pipeline and return structured run details."""
    validate_runtime_paths(args)
    input_path = Path(args.input)
    output_dir = Path(args.out_dir)
    input_variants = parse_vcf(input_path, args.assembly)

    clinvar_index = load_clinvar_index(
        variant_summary_path=Path(args.variant_summary),
        conflict_summary_path=Path(args.conflict_summary) if args.conflict_summary else None,
        submission_summary_path=Path(args.submission_summary) if args.submission_summary else None,
        target_variant_keys={variant.variant_key for variant in input_variants},
        cache_db_path=Path(args.clinvar_cache_db) if args.clinvar_cache_db else None,
        use_processed_cache=not bool(args.disable_clinvar_cache),
    )

    annotated_variants = annotate_variants(input_variants, clinvar_index)
    if args.enable_pharmgkb:
        pharmgkb_client = PharmGKBClient()
        annotated_variants = enrich_annotated_variants(annotated_variants, pharmgkb_client)
    else:
        pharmgkb_client = None
    ranked_variants = rank_variants(annotated_variants)

    run_metadata = build_run_metadata(args, output_dir, clinvar_index.provenance)
    if pharmgkb_client is not None:
        run_metadata.sources.extend(pharmgkb_client.provenance)
    run_metadata.statistics.input_variant_count = len(input_variants)
    run_metadata.statistics.clinvar_matched_count = sum(1 for item in annotated_variants if item.has_clinvar_match)
    run_metadata.statistics.clinvar_unmatched_count = len(annotated_variants) - run_metadata.statistics.clinvar_matched_count
    run_metadata.statistics.conflict_flagged_count = sum(1 for item in annotated_variants if item.has_conflict)
    run_metadata.statistics.pharmgkb_enriched_count = sum(
        1
        for item in annotated_variants
        if item.pharmgkb is not None and item.pharmgkb.matched
    )

    report_context = build_report_context(ranked_variants, run_metadata=run_metadata)
    variant_export_records = build_variant_export_records(ranked_variants)
    prioritized_variants_artifact = PrioritizedVariantsArtifact(records=variant_export_records)
    summary_artifact = SummaryArtifact(**report_context["summary_artifact"])
    csv_rows = [record.model_dump(mode="json") for record in variant_export_records]

    outputs = {
        "annotated_variants_csv": write_csv(output_dir / "annotated_variants.csv", csv_rows),
        "prioritized_variants_json": write_json(
            output_dir / "prioritized_variants.json",
            prioritized_variants_artifact.model_dump(mode="json"),
        ),
        "summary_json": write_json(output_dir / "summary.json", summary_artifact.model_dump(mode="json")),
        "run_metadata_json": write_json(output_dir / "run_metadata.json", run_metadata.model_dump(mode="json")),
        "report_html": write_html_report(output_dir / "report.html", ranked_variants, run_metadata=run_metadata),
    }
    return PipelineRunResult(outputs=outputs, run_metadata=run_metadata, report_context=report_context)


def run_pipeline_with_details(args: argparse.Namespace) -> tuple[dict[str, Path], RunMetadata]:
    """Execute the local ClinVar-first pipeline and return outputs with run metadata."""
    result = run_pipeline_with_result(args)
    return result.outputs, result.run_metadata


def run_pipeline(args: argparse.Namespace) -> dict[str, Path]:
    """Execute the local ClinVar-first pipeline and write output artifacts."""
    outputs, _ = run_pipeline_with_details(args)
    return outputs
