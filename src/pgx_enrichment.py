"""Optional pharmacogenomics enrichment helpers."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from time import sleep
from typing import Any

import requests

from .models import AnnotatedVariant, PharmGKBAnnotation

PHARMGKB_API_BASE = "https://api.pharmgkb.org/v1"
PHARMGKB_DOCS_URL = "https://api.pharmgkb.org/"
DEFAULT_CACHE_DIR = Path("data/pharmgkb/cache")
DEFAULT_TIMEOUT_SECONDS = 20
DEFAULT_RATE_LIMIT_DELAY_SECONDS = 0.55

PHARMGKB_ENDPOINTS = {
    "genes": f"{PHARMGKB_API_BASE}/data/gene",
    "chemicals": f"{PHARMGKB_API_BASE}/data/chemical",
    "variants": f"{PHARMGKB_API_BASE}/data/variant",
    "clinical_annotations": f"{PHARMGKB_API_BASE}/data/clinicalAnnotation",
    "guideline_annotations": f"{PHARMGKB_API_BASE}/data/guidelineAnnotation",
}


class PharmGKBClient:
    """Small cache-backed client for public PharmGKB enrichment queries."""

    def __init__(
        self,
        cache_dir: Path | None = None,
        timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
        session: requests.Session | None = None,
        rate_limit_delay_seconds: float = DEFAULT_RATE_LIMIT_DELAY_SECONDS,
    ) -> None:
        self.cache_dir = cache_dir or DEFAULT_CACHE_DIR
        self.timeout_seconds = timeout_seconds
        self.session = session or requests.Session()
        self.rate_limit_delay_seconds = rate_limit_delay_seconds
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _cache_path(self, endpoint_name: str, params: dict[str, str]) -> Path:
        """Build a deterministic cache filename for an endpoint + query."""
        serialized = json.dumps(
            {
                "endpoint_name": endpoint_name,
                "params": {key: params[key] for key in sorted(params)},
            },
            sort_keys=True,
        )
        digest = hashlib.sha256(serialized.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{endpoint_name}_{digest}.json"

    def _load_cache(self, cache_path: Path) -> dict[str, Any] | None:
        """Load a cached payload if present."""
        if not cache_path.exists():
            return None
        return json.loads(cache_path.read_text(encoding="utf-8"))

    def _write_cache(self, cache_path: Path, payload: dict[str, Any]) -> None:
        """Persist an API payload to the local cache."""
        cache_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    def _get(self, endpoint_name: str, params: dict[str, str]) -> tuple[list[dict[str, Any]], bool]:
        """Return a PharmGKB list payload, using cache and soft-fail behavior."""
        cache_path = self._cache_path(endpoint_name, params)
        cached = self._load_cache(cache_path)
        if cached is not None:
            return list(cached.get("data", [])), True

        endpoint = PHARMGKB_ENDPOINTS[endpoint_name]
        try:
            response = self.session.get(endpoint, params=params, timeout=self.timeout_seconds)
        except requests.RequestException:
            return [], False

        if response.status_code not in {200, 404}:
            return [], False

        try:
            payload = response.json()
        except ValueError:
            return [], False

        if response.status_code == 404 or payload.get("status") != "success":
            payload = {"status": "success", "data": []}

        self._write_cache(cache_path, payload)
        sleep(self.rate_limit_delay_seconds)
        return list(payload.get("data", [])), False

    def fetch_gene(self, symbol: str) -> tuple[list[dict[str, Any]], bool]:
        """Query PharmGKB genes by symbol."""
        return self._get("genes", {"symbol": symbol})

    def fetch_variant(self, symbol: str) -> tuple[list[dict[str, Any]], bool]:
        """Query PharmGKB variants by symbol."""
        return self._get("variants", {"symbol": symbol})

    def fetch_clinical_annotations_for_gene(self, symbol: str) -> tuple[list[dict[str, Any]], bool]:
        """Query PharmGKB clinical annotations by gene symbol."""
        return self._get("clinical_annotations", {"location.genes.symbol": symbol})

    def fetch_guideline_annotations_for_gene(self, symbol: str) -> tuple[list[dict[str, Any]], bool]:
        """Query PharmGKB guideline annotations by gene symbol."""
        return self._get("guideline_annotations", {"relatedGenes.symbol": symbol})


def _extract_gene_symbol(annotated_variant: AnnotatedVariant) -> str | None:
    """Choose the best available gene symbol for PharmGKB lookups."""
    if annotated_variant.clinvar.gene:
        return annotated_variant.clinvar.gene
    return annotated_variant.input_variant.gene


def enrich_annotated_variant(
    annotated_variant: AnnotatedVariant,
    client: PharmGKBClient,
) -> AnnotatedVariant:
    """Attach optional PharmGKB context to an annotated variant without raising on API failure."""
    gene_symbol = _extract_gene_symbol(annotated_variant)
    variant_symbol = annotated_variant.input_variant.variant_id

    from_cache = False
    gene_hits: list[dict[str, Any]] = []
    variant_hits: list[dict[str, Any]] = []
    clinical_hits: list[dict[str, Any]] = []
    guideline_hits: list[dict[str, Any]] = []

    if gene_symbol:
        gene_hits, gene_cached = client.fetch_gene(gene_symbol)
        clinical_hits, clinical_cached = client.fetch_clinical_annotations_for_gene(gene_symbol)
        guideline_hits, guideline_cached = client.fetch_guideline_annotations_for_gene(gene_symbol)
        from_cache = from_cache or gene_cached or clinical_cached or guideline_cached

    if variant_symbol:
        variant_hits, variant_cached = client.fetch_variant(variant_symbol)
        from_cache = from_cache or variant_cached

    annotation = PharmGKBAnnotation(
        queried=bool(gene_symbol or variant_symbol),
        matched=bool(gene_hits or variant_hits or clinical_hits or guideline_hits),
        from_cache=from_cache,
        gene_symbols=sorted({item.get("symbol", "") for item in gene_hits if item.get("symbol")}),
        pharmgkb_gene_ids=sorted({item.get("id", "") for item in gene_hits if item.get("id")}),
        pharmgkb_variant_ids=sorted({item.get("id", "") for item in variant_hits if item.get("id")}),
        chemicals=sorted(
            {
                chemical.get("name", "")
                for item in clinical_hits
                for chemical in item.get("relatedChemicals", [])
                if chemical.get("name")
            }
        ),
        clinical_annotation_ids=sorted(
            {item.get("accessionId", "") for item in clinical_hits if item.get("accessionId")}
        ),
        guideline_annotation_ids=sorted(
            {item.get("id", "") for item in guideline_hits if item.get("id")}
        ),
        evidence_notes=[],
    )

    if clinical_hits:
        annotation.evidence_notes.append(
            f"Found {len(clinical_hits)} PharmGKB clinical annotation record(s) for gene symbol {gene_symbol}."
        )
    if guideline_hits:
        annotation.evidence_notes.append(
            f"Found {len(guideline_hits)} PharmGKB guideline annotation record(s) for gene symbol {gene_symbol}."
        )
    if variant_hits and variant_symbol:
        annotation.evidence_notes.append(
            f"Found {len(variant_hits)} PharmGKB variant record(s) for symbol {variant_symbol}."
        )

    enriched = annotated_variant.model_copy(deep=True)
    enriched.pharmgkb = annotation
    if annotation.matched:
        enriched.flags.append("pharmgkb_enriched")
    elif annotation.queried:
        enriched.flags.append("pharmgkb_queried_no_match")
    return enriched


def enrich_annotated_variants(
    annotated_variants: list[AnnotatedVariant],
    client: PharmGKBClient,
) -> list[AnnotatedVariant]:
    """Enrich a batch of annotated variants with optional PharmGKB context."""
    return [enrich_annotated_variant(annotated_variant, client) for annotated_variant in annotated_variants]
