"""ClinVar snapshot loading and indexing helpers."""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, TextIO

import pandas as pd

from .models import (
    ClinVarMatch,
    ConflictSummary,
    DataProvenance,
    GenomeAssembly,
    InputVariant,
    MatchStrategy,
    SubmissionEvidence,
)
from .vcf_parser import normalize_allele, normalize_chromosome

VARIANT_SUMMARY_COLUMNS = [
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


def review_status_to_stars(review_status: str | None) -> int | None:
    """Convert ClinVar review status text into the familiar star count."""
    if not review_status:
        return None

    normalized = review_status.strip().lower()
    if normalized == "practice guideline":
        return 4
    if normalized == "reviewed by expert panel":
        return 3
    if normalized == "criteria provided, multiple submitters, no conflicts":
        return 2
    if normalized == "criteria provided, single submitter":
        return 1
    return 0


def _open_text_stream(path: Path) -> TextIO:
    """Open a plain-text or gzipped tabular ClinVar file."""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _parse_assembly(value: str | None) -> GenomeAssembly:
    """Normalize assembly labels from ClinVar files into supported enum values."""
    if value == GenomeAssembly.GRCH37.value:
        return GenomeAssembly.GRCH37
    if value == GenomeAssembly.GRCH38.value:
        return GenomeAssembly.GRCH38
    return GenomeAssembly.UNKNOWN


def _split_pipe_values(value: str | None) -> list[str]:
    """Split ClinVar multi-value fields while removing empty sentinels."""
    if value is None:
        return []
    results = []
    for part in value.split("|"):
        normalized = part.strip()
        if normalized and normalized.lower() not in {"-", "na", "not provided"}:
            results.append(normalized)
    return results


def _choose_preferred_match(existing: ClinVarMatch, candidate: ClinVarMatch) -> ClinVarMatch:
    """Prefer the stronger ClinVar aggregate when duplicate coordinate keys appear."""
    existing_score = (
        existing.review_stars or -1,
        1 if existing.clinical_significance else 0,
        existing.variation_id or -1,
    )
    candidate_score = (
        candidate.review_stars or -1,
        1 if candidate.clinical_significance else 0,
        candidate.variation_id or -1,
    )
    return candidate if candidate_score > existing_score else existing


def _iter_submission_rows(submission_summary_path: Path) -> Iterator[dict[str, str]]:
    """Yield parsed submission-summary rows after skipping the explanatory preamble."""
    with _open_text_stream(submission_summary_path) as handle:
        header_line = None
        for line in handle:
            if line.startswith("#VariationID\t"):
                header_line = line.lstrip("#").rstrip("\n")
                break

        if header_line is None:
            raise ValueError("Submission summary file does not contain a tabular header.")

        fieldnames = header_line.split("\t")
        reader = csv.DictReader(handle, fieldnames=fieldnames, delimiter="\t")
        for row in reader:
            yield {key: (value or "").strip() for key, value in row.items()}


def _build_provenance(source_name: str, source_kind: str, source_path: Path) -> DataProvenance:
    """Construct provenance metadata from a local file path."""
    return DataProvenance(
        source_name=source_name,
        source_kind=source_kind,
        source_path=str(source_path),
    )


def _parse_int_field(value: str | None) -> int | None:
    """Parse an integer field from ClinVar text, returning None for blanks."""
    if value is None:
        return None
    normalized = value.strip()
    if not normalized:
        return None
    try:
        return int(normalized)
    except ValueError:
        return None


@dataclass(slots=True)
class ClinVarIndex:
    """In-memory ClinVar lookup tables used by the annotation stage."""

    exact_matches: dict[tuple[str, str, int, str, str], ClinVarMatch]
    conflicts_by_variation_id: dict[int, ConflictSummary] = field(default_factory=dict)
    submissions_by_variation_id: dict[int, SubmissionEvidence] = field(default_factory=dict)
    provenance: list[DataProvenance] = field(default_factory=list)

    def lookup(self, input_variant: InputVariant) -> ClinVarMatch:
        """Return the best exact ClinVar match for a normalized input variant."""
        match = self.exact_matches.get(input_variant.variant_key)
        if match is None:
            return ClinVarMatch()

        resolved = match.model_copy(deep=True)
        variation_id = resolved.variation_id
        if variation_id is not None:
            conflict = self.conflicts_by_variation_id.get(variation_id)
            if conflict is not None:
                resolved.conflict = conflict.model_copy(deep=True)
            submission = self.submissions_by_variation_id.get(variation_id)
            if submission is not None:
                resolved.submissions = submission.model_copy(deep=True)
        return resolved


def load_variant_summary_index(
    variant_summary_path: Path,
    chunk_size: int = 50_000,
) -> ClinVarIndex:
    """Stream `variant_summary.txt.gz` into an exact coordinate lookup index."""
    exact_matches: dict[tuple[str, str, int, str, str], ClinVarMatch] = {}

    reader = pd.read_csv(
        variant_summary_path,
        sep="\t",
        compression="gzip" if variant_summary_path.suffix.lower() == ".gz" else None,
        usecols=VARIANT_SUMMARY_COLUMNS,
        dtype=str,
        chunksize=chunk_size,
        low_memory=False,
    )

    for chunk in reader:
        chunk = chunk.fillna("")
        for row in chunk.itertuples(index=False):
            assembly = _parse_assembly(row.Assembly)
            if assembly == GenomeAssembly.UNKNOWN:
                continue

            position_vcf = (row.PositionVCF or "").strip()
            reference_allele = normalize_allele((row.ReferenceAlleleVCF or "").strip())
            alternate_allele = normalize_allele((row.AlternateAlleleVCF or "").strip())
            chromosome = normalize_chromosome((row.Chromosome or "").strip())
            if not position_vcf or not chromosome:
                continue
            if reference_allele in {"", "-"} or alternate_allele in {"", "-"}:
                continue

            variation_id = _parse_int_field(row.VariationID)
            allele_id = _parse_int_field(row[0])
            parsed_position = _parse_int_field(position_vcf)
            if variation_id is None or allele_id is None or parsed_position is None:
                continue
            candidate = ClinVarMatch(
                matched=True,
                match_strategy=MatchStrategy.EXACT,
                assembly=assembly,
                chromosome=chromosome,
                position=parsed_position,
                reference_allele=reference_allele,
                alternate_allele=alternate_allele,
                variation_id=variation_id,
                allele_id=allele_id,
                accession=(row.RCVaccession or "").split("|", 1)[0] or None,
                preferred_name=(row.Name or "").strip() or None,
                gene=(row.GeneSymbol or "").strip() or None,
                condition_names=_split_pipe_values((row.PhenotypeList or "").strip()),
                clinical_significance=(row.ClinicalSignificance or "").strip() or None,
                review_status=(row.ReviewStatus or "").strip() or None,
                review_stars=review_status_to_stars((row.ReviewStatus or "").strip() or None),
                interpretation_origin=(row.Origin or "").strip() or None,
                last_evaluated=(row.LastEvaluated or "").strip() or None,
            )
            key = candidate.variant_key
            existing = exact_matches.get(key)
            exact_matches[key] = candidate if existing is None else _choose_preferred_match(existing, candidate)

    return ClinVarIndex(
        exact_matches=exact_matches,
        provenance=[_build_provenance("ClinVar variant summary", "file", variant_summary_path)],
    )


def load_conflict_lookup(
    conflict_summary_path: Path,
    target_variation_ids: set[int] | None = None,
) -> tuple[dict[int, ConflictSummary], DataProvenance]:
    """Load conflict summaries keyed by ClinVar VariationID."""
    accumulators: dict[int, dict[str, set[str]]] = {}

    with _open_text_stream(conflict_summary_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            variation_text = (row.get("NCBI_Variation_ID") or "").strip()
            if not variation_text:
                continue

            variation_id = int(variation_text)
            if target_variation_ids is not None and variation_id not in target_variation_ids:
                continue

            bucket = accumulators.setdefault(
                variation_id,
                {"submitters": set(), "significances": set(), "preferred_names": set()},
            )
            for key in ("Submitter1", "Submitter2"):
                value = (row.get(key) or "").strip()
                if value:
                    bucket["submitters"].add(value)
            for key in ("Submitter1_ClinSig", "Submitter2_ClinSig"):
                value = (row.get(key) or "").strip()
                if value:
                    bucket["significances"].add(value)
            preferred = (row.get("ClinVar_Preferred") or "").strip()
            if preferred:
                bucket["preferred_names"].add(preferred)

    conflicts = {}
    for variation_id, bucket in accumulators.items():
        summary_text = None
        if bucket["preferred_names"]:
            summary_text = "; ".join(sorted(bucket["preferred_names"]))
        conflicts[variation_id] = ConflictSummary(
            has_conflict=True,
            conflict_significance=sorted(bucket["significances"]),
            submitter_count=len(bucket["submitters"]) or None,
            summary_text=summary_text,
        )
    return conflicts, _build_provenance("ClinVar conflicting interpretations", "file", conflict_summary_path)


def load_submission_lookup(
    submission_summary_path: Path,
    target_variation_ids: set[int] | None = None,
) -> tuple[dict[int, SubmissionEvidence], DataProvenance]:
    """Load submission aggregates keyed by ClinVar VariationID."""
    accumulators: dict[int, dict[str, int | set[str]]] = {}
    for row in _iter_submission_rows(submission_summary_path):
        variation_text = (row.get("VariationID") or "").strip()
        if not variation_text:
            continue

        variation_id = int(variation_text)
        if target_variation_ids is not None and variation_id not in target_variation_ids:
            continue

        bucket = accumulators.setdefault(
            variation_id,
            {
                "total_submissions": 0,
                "submitter_names": set(),
                "review_statuses": set(),
                "clinical_significances": set(),
            },
        )
        bucket["total_submissions"] = int(bucket["total_submissions"]) + 1

        submitter = (row.get("Submitter") or "").strip()
        if submitter:
            bucket["submitter_names"].add(submitter)

        review_status = (row.get("ReviewStatus") or "").strip()
        if review_status:
            bucket["review_statuses"].add(review_status)

        clinical_significance = (row.get("ClinicalSignificance") or "").strip()
        if clinical_significance:
            bucket["clinical_significances"].add(clinical_significance)

    submissions = {}
    for variation_id, bucket in accumulators.items():
        submissions[variation_id] = SubmissionEvidence(
            total_submissions=int(bucket["total_submissions"]),
            submitter_names=sorted(bucket["submitter_names"]),
            review_statuses=sorted(bucket["review_statuses"]),
            clinical_significances=sorted(bucket["clinical_significances"]),
        )
    return submissions, _build_provenance("ClinVar submission summary", "file", submission_summary_path)


def enrich_index_with_supporting_data(
    index: ClinVarIndex,
    conflict_summary_path: Path | None = None,
    submission_summary_path: Path | None = None,
    target_variation_ids: Iterable[int] | None = None,
) -> ClinVarIndex:
    """Attach conflict and submission layers to an existing exact-match index."""
    allowed_variation_ids = set(target_variation_ids) if target_variation_ids is not None else {
        match.variation_id for match in index.exact_matches.values() if match.variation_id is not None
    }

    if conflict_summary_path is not None:
        conflicts, provenance = load_conflict_lookup(conflict_summary_path, allowed_variation_ids)
        index.conflicts_by_variation_id.update(conflicts)
        index.provenance.append(provenance)

    if submission_summary_path is not None:
        submissions, provenance = load_submission_lookup(submission_summary_path, allowed_variation_ids)
        index.submissions_by_variation_id.update(submissions)
        index.provenance.append(provenance)

    return index


def load_clinvar_index(
    variant_summary_path: Path,
    conflict_summary_path: Path | None = None,
    submission_summary_path: Path | None = None,
    target_variation_ids: Iterable[int] | None = None,
    chunk_size: int = 50_000,
) -> ClinVarIndex:
    """Build the exact-match index and optionally attach supporting evidence layers."""
    index = load_variant_summary_index(variant_summary_path, chunk_size=chunk_size)
    return enrich_index_with_supporting_data(
        index,
        conflict_summary_path=conflict_summary_path,
        submission_summary_path=submission_summary_path,
        target_variation_ids=target_variation_ids,
    )
