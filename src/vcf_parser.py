"""VCF parsing utilities."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import TextIO

from .models import GenomeAssembly, InputVariant

VCF_REQUIRED_COLUMNS = (
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
)

GENE_INFO_KEYS = ("GENE", "GENE_SYMBOL", "SYMBOL", "HGNC")
TRANSCRIPT_INFO_KEYS = ("TRANSCRIPT", "FEATURE", "FEATUREID")
CONSEQUENCE_INFO_KEYS = ("CONSEQUENCE", "CSQ_CONSEQUENCE")
IMPACT_INFO_KEYS = ("IMPACT",)
MISSING_FIELD_VALUES = {"", ".", "-", "na", "not provided"}


def open_variant_stream(input_path: Path) -> TextIO:
    """Open a plain-text or gzipped VCF input file."""
    if input_path.suffix.lower() == ".gz":
        return gzip.open(input_path, "rt", encoding="utf-8")
    return input_path.open("r", encoding="utf-8")


def normalize_chromosome(chromosome: str) -> str:
    """Normalize common chromosome labels to a ClinVar-friendly form."""
    normalized = chromosome.strip()
    if normalized.lower().startswith("chr"):
        normalized = normalized[3:]
    if normalized.upper() == "M":
        return "MT"
    return normalized


def normalize_allele(allele: str) -> str:
    """Normalize allele text for deterministic coordinate matching."""
    return allele.strip().upper()


def parse_info_field(info_field: str) -> dict[str, str]:
    """Parse key-value entries from a VCF INFO field."""
    if not info_field or info_field == ".":
        return {}

    parsed: dict[str, str] = {}
    for entry in info_field.split(";"):
        if not entry:
            continue
        if "=" not in entry:
            parsed[entry] = "true"
            continue
        key, value = entry.split("=", 1)
        parsed[key] = value
    return parsed


def _extract_gene(info_map: dict[str, str]) -> str | None:
    """Extract a gene symbol from direct INFO keys or common annotation formats."""
    for key in GENE_INFO_KEYS:
        value = info_map.get(key)
        if value and value.strip().lower() not in MISSING_FIELD_VALUES:
            return value.split("|", 1)[0]

    gene_info = info_map.get("GENEINFO")
    if gene_info and gene_info.strip().lower() not in MISSING_FIELD_VALUES:
        return gene_info.split(":", 1)[0]

    ann_fields = _parse_ann_fields(info_map)
    if ann_fields.get("gene"):
        return ann_fields["gene"]
    return None


def _extract_transcript(info_map: dict[str, str]) -> str | None:
    """Extract a transcript identifier from common INFO keys or SnpEff ANN."""
    for key in TRANSCRIPT_INFO_KEYS:
        value = info_map.get(key)
        if value and value.strip().lower() not in MISSING_FIELD_VALUES:
            return value.split("|", 1)[0]

    ann_fields = _parse_ann_fields(info_map)
    if ann_fields.get("transcript"):
        return ann_fields["transcript"]
    return None


def _extract_consequence(info_map: dict[str, str]) -> str | None:
    """Extract consequence text from direct keys or SnpEff ANN."""
    for key in CONSEQUENCE_INFO_KEYS:
        value = info_map.get(key)
        if value and value.strip().lower() not in MISSING_FIELD_VALUES:
            return value

    ann_fields = _parse_ann_fields(info_map)
    if ann_fields.get("consequence"):
        return ann_fields["consequence"]
    return None


def _extract_impact(info_map: dict[str, str]) -> str | None:
    """Extract impact text from direct keys or SnpEff ANN."""
    for key in IMPACT_INFO_KEYS:
        value = info_map.get(key)
        if value and value.strip().lower() not in MISSING_FIELD_VALUES:
            return value

    ann_fields = _parse_ann_fields(info_map)
    if ann_fields.get("impact"):
        return ann_fields["impact"]
    return None


def _parse_ann_fields(info_map: dict[str, str]) -> dict[str, str]:
    """Extract selected SnpEff ANN fields from the first annotation entry."""
    ann_value = info_map.get("ANN")
    if not ann_value:
        return {}

    first_annotation = ann_value.split(",", 1)[0]
    fields = first_annotation.split("|")
    if len(fields) < 7:
        return {}

    return {
        "consequence": fields[1].strip() or "",
        "impact": fields[2].strip() or "",
        "gene": fields[3].strip() or "",
        "transcript": fields[6].strip() or "",
    }


def _split_vcf_fields(line: str) -> list[str]:
    """Split a VCF line, preferring tab delimiters but tolerating whitespace."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) >= 8:
        return fields
    return line.split()


def _validate_header(header_fields: list[str]) -> None:
    """Ensure the VCF header exposes the required leading columns."""
    if len(header_fields) < len(VCF_REQUIRED_COLUMNS):
        raise ValueError("VCF header has fewer than 8 required columns.")

    leading_columns = tuple(header_fields[: len(VCF_REQUIRED_COLUMNS)])
    if leading_columns != VCF_REQUIRED_COLUMNS:
        raise ValueError(
            "VCF header does not match required leading columns: "
            f"{', '.join(VCF_REQUIRED_COLUMNS)}"
        )


def _build_input_variant(
    record_index: int,
    assembly: GenomeAssembly,
    chrom: str,
    pos: str,
    variant_id: str,
    ref: str,
    alt: str,
    qual: str,
    filter_status: str,
    info: str,
) -> InputVariant:
    """Convert parsed VCF fields into a normalized InputVariant."""
    info_map = parse_info_field(info)
    normalized_id = None if variant_id == "." else variant_id

    return InputVariant(
        record_id=f"record-{record_index}",
        assembly=assembly,
        chromosome=normalize_chromosome(chrom),
        position=int(pos),
        reference_allele=normalize_allele(ref),
        alternate_allele=normalize_allele(alt),
        quality=None if qual == "." else qual,
        filter_status=None if filter_status == "." else filter_status,
        variant_id=normalized_id,
        gene=_extract_gene(info_map),
        transcript=_extract_transcript(info_map),
        consequence=_extract_consequence(info_map),
        impact=_extract_impact(info_map),
        info=info_map,
    )


def parse_vcf(input_path: Path, assembly: GenomeAssembly, max_variants: int | None = None) -> list[InputVariant]:
    """Read a VCF or VCF.GZ file and return normalized InputVariant records."""
    saw_header = False
    variants: list[InputVariant] = []

    with open_variant_stream(input_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header_fields = _split_vcf_fields(line)
                _validate_header(header_fields)
                saw_header = True
                continue

            if line.startswith("#"):
                continue

            if not saw_header:
                raise ValueError("VCF record encountered before #CHROM header line.")

            fields = _split_vcf_fields(line)
            if len(fields) < 8:
                raise ValueError(f"VCF record has fewer than 8 columns at line {line_number}.")

            chrom, pos, variant_id, ref, alt_field, qual, filter_status, info = fields[:8]
            for alt_allele in alt_field.split(","):
                variants.append(
                    _build_input_variant(
                        record_index=len(variants) + 1,
                        assembly=assembly,
                        chrom=chrom,
                        pos=pos,
                        variant_id=variant_id,
                        ref=ref,
                        alt=alt_allele,
                        qual=qual,
                        filter_status=filter_status,
                        info=info,
                    )
                )
                if max_variants is not None and len(variants) >= max_variants:
                    return variants

    if not saw_header:
        raise ValueError("Input does not contain a #CHROM header line.")

    return variants
