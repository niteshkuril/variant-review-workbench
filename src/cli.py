"""CLI entry point for variant-review-workbench."""

from __future__ import annotations

import argparse
from pathlib import Path

from .app_service import PipelineUsageError, run_pipeline, run_pipeline_with_details
from .models import GenomeAssembly, RunMetadata


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


def main() -> None:
    """Parse CLI arguments and execute the local reporting pipeline."""
    parser = build_parser()
    args = parser.parse_args()
    try:
        outputs, run_metadata = run_pipeline_with_details(args)
    except PipelineUsageError as error:
        parser.exit(2, f"Error: {error}\n")
    except ValueError as error:
        parser.exit(2, f"Error: {error}\n")
    except OSError as error:
        parser.exit(1, f"Runtime error: {error}\n")

    _emit_completion_summary(outputs, run_metadata)


if __name__ == "__main__":
    main()
