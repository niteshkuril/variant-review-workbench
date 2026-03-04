"""Variant annotation orchestration."""

from __future__ import annotations

from .clinvar_index import ClinVarIndex
from .models import AnnotatedVariant, InputVariant


def _normalize_gene_symbol(symbol: str | None) -> str | None:
    """Normalize gene symbols before comparing user and ClinVar annotations."""
    if symbol is None:
        return None
    normalized = symbol.strip().upper()
    return normalized or None


def build_annotation_flags(input_variant: InputVariant, annotated_variant: AnnotatedVariant) -> list[str]:
    """Derive transparent workflow flags from the merged annotation state."""
    flags: list[str] = []
    clinvar = annotated_variant.clinvar

    if not clinvar.matched:
        flags.append("clinvar_unmatched")
        return flags

    flags.append("clinvar_matched")

    if clinvar.conflict.has_conflict:
        flags.append("clinvar_conflict")

    if clinvar.review_stars is not None:
        flags.append(f"clinvar_review_stars_{clinvar.review_stars}")

    input_gene = _normalize_gene_symbol(input_variant.gene)
    clinvar_gene = _normalize_gene_symbol(clinvar.gene)
    if input_gene and clinvar_gene and input_gene != clinvar_gene:
        flags.append("gene_symbol_mismatch")

    if clinvar.submissions is not None and clinvar.submissions.total_submissions:
        flags.append("submission_evidence_available")

    return flags


def annotate_variant(input_variant: InputVariant, clinvar_index: ClinVarIndex) -> AnnotatedVariant:
    """Attach ClinVar context to a single normalized input variant."""
    clinvar_match = clinvar_index.lookup(input_variant)
    annotated = AnnotatedVariant(
        input_variant=input_variant,
        clinvar=clinvar_match,
    )
    annotated.flags = build_annotation_flags(input_variant, annotated)
    return annotated


def annotate_variants(
    input_variants: list[InputVariant],
    clinvar_index: ClinVarIndex,
) -> list[AnnotatedVariant]:
    """Annotate a batch of normalized variants using the local ClinVar index."""
    return [annotate_variant(input_variant, clinvar_index) for input_variant in input_variants]
