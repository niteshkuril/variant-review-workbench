"""CLI entry point for variant-review-workbench."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .annotator import annotate_variants
from .clinvar_index import load_clinvar_index
from .models import (
    GenomeAssembly,
    PrioritizedVariantsArtifact,
    RunMetadata,
    SummaryArtifact,
    VariantExportRecord,
)
from .pgx_enrichment import PharmGKBClient, enrich_annotated_variants
from .ranker import rank_variants
from .report_builder import build_report_context, write_html_report
from .vcf_parser import parse_vcf


class CliUsageError(ValueError):
    """Raised when CLI arguments are syntactically valid but operationally unusable."""


def _parse_assembly(value: str) -> GenomeAssembly:
    """Parse a CLI assembly argument into a supported enum value."""
    normalized = value.strip()
    if normalized == GenomeAssembly.GRCH37.value:
        return GenomeAssembly.GRCH37
    if normalized == GenomeAssembly.GRCH38.value:
        return GenomeAssembly.GRCH38
    raise argparse.ArgumentTypeError("assembly must be either 'GRCh37' or 'GRCh38'")


def build_parser() -> argparse.ArgumentParser:
    """Create the top-level command-line parser."""
    parser = argparse.ArgumentParser(
        description="Annotate, rank, and report small variants against a local ClinVar snapshot.",
    )
    parser.add_argument("--input", required=True, help="Path to the input VCF or VCF.GZ file.")
    parser.add_argument("--assembly", required=True, type=_parse_assembly, help="Reference assembly: GRCh37 or GRCh38.")
    parser.add_argument("--variant-summary", required=True, help="Path to ClinVar variant_summary.txt.gz.")
    parser.add_argument(
        "--conflict-summary",
        default=None,
        help="Optional path to summary_of_conflicting_interpretations.txt.",
    )
    parser.add_argument(
        "--submission-summary",
        default=None,
        help="Optional path to submission_summary.txt.gz.",
    )
    parser.add_argument(
        "--clinvar-cache-db",
        default=None,
        help="Optional path to a persistent processed ClinVar SQLite cache.",
    )
    parser.add_argument(
        "--disable-clinvar-cache",
        action="store_true",
        help="Disable the persistent processed ClinVar cache and read raw files directly.",
    )
    parser.add_argument("--out-dir", required=True, help="Directory where run outputs will be written.")
    parser.add_argument(
        "--enable-pharmgkb",
        action="store_true",
        help="Enable optional live PharmGKB enrichment with local response caching.",
    )
    return parser


def _build_run_metadata(args: argparse.Namespace, output_dir: Path, index_sources: list) -> RunMetadata:
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
        raise CliUsageError(f"{label} was not found: {path}")
    if not path.is_file():
        raise CliUsageError(f"{label} must be a file path: {path}")


def _validate_runtime_paths(args: argparse.Namespace) -> None:
    """Validate user-supplied paths before executing the pipeline."""
    _validate_existing_file(Path(args.input), "Input VCF")
    _validate_existing_file(Path(args.variant_summary), "ClinVar variant summary")
    if args.conflict_summary:
        _validate_existing_file(Path(args.conflict_summary), "ClinVar conflict summary")
    if args.submission_summary:
        _validate_existing_file(Path(args.submission_summary), "ClinVar submission summary")

    output_dir = Path(args.out_dir)
    if output_dir.exists() and not output_dir.is_dir():
        raise CliUsageError(f"Output path must be a directory: {output_dir}")

    if args.clinvar_cache_db:
        cache_db_path = Path(args.clinvar_cache_db)
        if cache_db_path.exists() and cache_db_path.is_dir():
            raise CliUsageError(f"ClinVar cache path must be a file path, not a directory: {cache_db_path}")


def _emit_completion_summary(
    outputs: dict[str, Path],
    run_metadata: RunMetadata,
) -> None:
    """Print concise completion details and output locations for successful runs."""
    stats = run_metadata.statistics
    print(
        (
            "Run completed: "
            f"{stats.input_variant_count} input variant(s), "
            f"{stats.clinvar_matched_count} ClinVar match(es), "
            f"{stats.conflict_flagged_count} conflict-flagged, "
            f"{stats.pharmgkb_enriched_count} PharmGKB-enriched."
        )
    )
    if run_metadata.pharmgkb_enabled and stats.pharmgkb_enriched_count == 0:
        print("PharmGKB was enabled but no enrichment matches were found.")
    for label, path in outputs.items():
        print(f"{label}: {path.resolve()}")


def _build_variant_export_records(ranked_variants: list) -> list[VariantExportRecord]:
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


def _write_json(output_path: Path, payload: object) -> Path:
    """Write a JSON artifact to disk."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return output_path


def _write_csv(output_path: Path, rows: list[dict[str, object]]) -> Path:
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


def run_pipeline_with_details(args: argparse.Namespace) -> tuple[dict[str, Path], RunMetadata]:
    """Execute the local ClinVar-first pipeline and return outputs with run metadata."""
    _validate_runtime_paths(args)
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

    run_metadata = _build_run_metadata(args, output_dir, clinvar_index.provenance)
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
    variant_export_records = _build_variant_export_records(ranked_variants)
    prioritized_variants_artifact = PrioritizedVariantsArtifact(records=variant_export_records)
    summary_artifact = SummaryArtifact(**report_context["summary_artifact"])
    csv_rows = [record.model_dump(mode="json") for record in variant_export_records]

    outputs = {
        "annotated_variants_csv": _write_csv(output_dir / "annotated_variants.csv", csv_rows),
        "prioritized_variants_json": _write_json(
            output_dir / "prioritized_variants.json",
            prioritized_variants_artifact.model_dump(mode="json"),
        ),
        "summary_json": _write_json(output_dir / "summary.json", summary_artifact.model_dump(mode="json")),
        "run_metadata_json": _write_json(output_dir / "run_metadata.json", run_metadata.model_dump(mode="json")),
        "report_html": write_html_report(output_dir / "report.html", ranked_variants, run_metadata=run_metadata),
    }
    return outputs, run_metadata


def run_pipeline(args: argparse.Namespace) -> dict[str, Path]:
    """Execute the local ClinVar-first pipeline and write output artifacts."""
    outputs, _ = run_pipeline_with_details(args)
    return outputs


def main() -> None:
    """Parse CLI arguments and execute the local reporting pipeline."""
    parser = build_parser()
    args = parser.parse_args()
    try:
        outputs, run_metadata = run_pipeline_with_details(args)
    except CliUsageError as error:
        parser.exit(2, f"Error: {error}\n")
    except ValueError as error:
        parser.exit(2, f"Error: {error}\n")
    except OSError as error:
        parser.exit(1, f"Runtime error: {error}\n")

    _emit_completion_summary(outputs, run_metadata)


if __name__ == "__main__":
    main()
