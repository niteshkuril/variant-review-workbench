"""CLI entry point for variant-review-workbench."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .annotator import annotate_variants
from .clinvar_index import load_clinvar_index
from .models import GenomeAssembly, RunMetadata
from .ranker import rank_variants
from .report_builder import build_report_context, write_html_report
from .vcf_parser import parse_vcf


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
        description="Annotate and rank small variants against a local ClinVar snapshot.",
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
    parser.add_argument("--out-dir", required=True, help="Directory where run outputs will be written.")
    parser.add_argument(
        "--enable-pharmgkb",
        action="store_true",
        help="Reserved flag for future PharmGKB enrichment integration.",
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


def _build_variant_export_rows(report_context: dict[str, object]) -> list[dict[str, object]]:
    """Extract report table rows for machine-readable output writing."""
    rows = []
    for row in report_context["variant_rows"]:
        rows.append(
            {
                "record_id": row["record_id"],
                "gene": row["gene"],
                "locus": row["locus"],
                "priority_tier": row["priority_tier"],
                "priority_score": row["priority_score"],
                "clinical_significance": row["clinical_significance"],
                "review_status": row["review_status"],
                "condition_names": row["condition_names"],
                "conflict": row["conflict"],
                "pharmgkb": row["pharmgkb"],
                "flags": list(row["flags"]),
                "rationale": list(row["rationale"]),
            }
        )
    return rows


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
    lines = [",".join(header)]
    for row in rows:
        serialized = []
        for key in header:
            value = row[key]
            if isinstance(value, list):
                cell = "|".join(str(item) for item in value)
            else:
                cell = str(value)
            if any(token in cell for token in [",", '"', "\n"]):
                cell = '"' + cell.replace('"', '""') + '"'
            serialized.append(cell)
        lines.append(",".join(serialized))

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def run_pipeline(args: argparse.Namespace) -> dict[str, Path]:
    """Execute the local ClinVar-first pipeline and write output artifacts."""
    input_path = Path(args.input)
    output_dir = Path(args.out_dir)

    clinvar_index = load_clinvar_index(
        variant_summary_path=Path(args.variant_summary),
        conflict_summary_path=Path(args.conflict_summary) if args.conflict_summary else None,
        submission_summary_path=Path(args.submission_summary) if args.submission_summary else None,
    )

    input_variants = parse_vcf(input_path, args.assembly)
    annotated_variants = annotate_variants(input_variants, clinvar_index)
    ranked_variants = rank_variants(annotated_variants)

    run_metadata = _build_run_metadata(args, output_dir, clinvar_index.provenance)
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
    variant_rows = _build_variant_export_rows(report_context)

    outputs = {
        "annotated_variants_csv": _write_csv(output_dir / "annotated_variants.csv", variant_rows),
        "prioritized_variants_json": _write_json(output_dir / "prioritized_variants.json", variant_rows),
        "summary_json": _write_json(output_dir / "summary.json", report_context["summary"]),
        "run_metadata_json": _write_json(output_dir / "run_metadata.json", run_metadata.model_dump(mode="json")),
        "report_html": write_html_report(output_dir / "report.html", ranked_variants, run_metadata=run_metadata),
    }
    return outputs


def main() -> None:
    """Parse CLI arguments and execute the local reporting pipeline."""
    args = build_parser().parse_args()
    outputs = run_pipeline(args)
    for label, path in outputs.items():
        print(f"{label}: {path.resolve()}")


if __name__ == "__main__":
    main()
