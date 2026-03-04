from __future__ import annotations

import argparse
import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

import requests

from src.cli import run_pipeline
from src.models import (
    AnnotatedVariant,
    ClinVarMatch,
    GenomeAssembly,
    InputVariant,
)
from src.pgx_enrichment import PharmGKBClient, enrich_annotated_variant


def build_annotated_variant() -> AnnotatedVariant:
    return AnnotatedVariant(
        input_variant=InputVariant(
            record_id="record-1",
            assembly=GenomeAssembly.GRCH38,
            chromosome="16",
            position=31105878,
            reference_allele="C",
            alternate_allele="T",
            variant_id="rs104894541",
            gene="VKORC1",
        ),
        clinvar=ClinVarMatch(
            matched=True,
            assembly=GenomeAssembly.GRCH38,
            chromosome="16",
            position=31105878,
            reference_allele="C",
            alternate_allele="T",
            gene="VKORC1",
            variation_id=1234,
        ),
        flags=["clinvar_matched"],
    )


class PharmGKBEnrichmentTests(unittest.TestCase):
    def test_client_caches_successful_responses(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            session = Mock()
            response = Mock()
            response.status_code = 200
            response.json.return_value = {
                "status": "success",
                "data": [{"objCls": "Gene", "id": "PA1", "symbol": "VKORC1"}],
            }
            session.get.return_value = response

            client = PharmGKBClient(
                cache_dir=Path(tmpdir),
                session=session,
                rate_limit_delay_seconds=0.0,
            )
            first_data, first_cached = client.fetch_gene("VKORC1")
            second_data, second_cached = client.fetch_gene("VKORC1")

        self.assertFalse(first_cached)
        self.assertTrue(second_cached)
        self.assertEqual(first_data, second_data)
        self.assertEqual(session.get.call_count, 1)

    def test_client_handles_request_failure_without_raising(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            session = Mock()
            session.get.side_effect = requests.RequestException("network down")
            client = PharmGKBClient(
                cache_dir=Path(tmpdir),
                session=session,
                rate_limit_delay_seconds=0.0,
            )

            data, cached = client.fetch_variant("rs104894541")

        self.assertEqual(data, [])
        self.assertFalse(cached)

    def test_enrich_annotated_variant_populates_pharmgkb_annotation(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            session = Mock()

            def fake_get(url: str, params: dict[str, str], timeout: int):
                response = Mock()
                response.status_code = 200
                if "gene" in url:
                    response.json.return_value = {"status": "success", "data": [{"id": "PA133787052", "symbol": "VKORC1"}]}
                elif "variant" in url:
                    response.json.return_value = {"status": "success", "data": [{"id": "PA166155139", "symbol": "rs104894541"}]}
                elif "clinicalAnnotation" in url:
                    response.json.return_value = {
                        "status": "success",
                        "data": [
                            {
                                "accessionId": "PA166134613",
                                "relatedChemicals": [{"name": "warfarin"}],
                            }
                        ],
                    }
                else:
                    response.json.return_value = {"status": "success", "data": [{"id": "PA166104949"}]}
                return response

            session.get.side_effect = fake_get
            client = PharmGKBClient(
                cache_dir=Path(tmpdir),
                session=session,
                rate_limit_delay_seconds=0.0,
            )

            enriched = enrich_annotated_variant(build_annotated_variant(), client)

        self.assertIsNotNone(enriched.pharmgkb)
        self.assertTrue(enriched.pharmgkb.matched)
        self.assertIn("VKORC1", enriched.pharmgkb.gene_symbols)
        self.assertIn("PA166155139", enriched.pharmgkb.pharmgkb_variant_ids)
        self.assertIn("warfarin", enriched.pharmgkb.chemicals)
        self.assertIn("pharmgkb_enriched", enriched.flags)

    def test_run_pipeline_with_pharmgkb_enabled_updates_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input.vcf"
            input_path.write_text(
                (
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "16\t31105878\trs104894541\tC\tT\t100\tPASS\tGENE=VKORC1;IMPACT=HIGH\n"
                ),
                encoding="utf-8",
            )

            variant_summary = root / "variant_summary.txt.gz"
            with gzip.open(variant_summary, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(
                    [
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
                )
                writer.writerow(
                    [
                        "10",
                        "single nucleotide variant",
                        "VKORC1 example",
                        "VKORC1",
                        "drug response",
                        "Jan 01, 2025",
                        "Warfarin response",
                        "criteria provided, single submitter",
                        "germline",
                        "GRCh38",
                        "16",
                        "1234",
                        "31105878",
                        "C",
                        "T",
                        "RCV000000001",
                    ]
                )

            fake_payloads = {
                "genes": [{"id": "PA133787052", "symbol": "VKORC1"}],
                "variants": [{"id": "PA166155139", "symbol": "rs104894541"}],
                "clinical_annotations": [{"accessionId": "PA166134613", "relatedChemicals": [{"name": "warfarin"}]}],
                "guideline_annotations": [{"id": "PA166104949"}],
            }

            original_get = PharmGKBClient._get

            def fake_get(self, endpoint_name: str, params: dict[str, str]):
                return fake_payloads.get(endpoint_name, []), False

            PharmGKBClient._get = fake_get
            try:
                out_dir = root / "outputs"
                args = argparse.Namespace(
                    input=str(input_path),
                    assembly=GenomeAssembly.GRCH38,
                    variant_summary=str(variant_summary),
                    conflict_summary=None,
                    submission_summary=None,
                    out_dir=str(out_dir),
                    enable_pharmgkb=True,
                )

                run_pipeline(args)
            finally:
                PharmGKBClient._get = original_get

            metadata = json.loads((out_dir / "run_metadata.json").read_text(encoding="utf-8"))
            prioritized = json.loads((out_dir / "prioritized_variants.json").read_text(encoding="utf-8"))

        self.assertTrue(metadata["pharmgkb_enabled"])
        self.assertEqual(metadata["statistics"]["pharmgkb_enriched_count"], 1)
        self.assertEqual(prioritized[0]["pharmgkb"], "Yes")
        self.assertIn("pharmgkb_enriched", prioritized[0]["flags"])


if __name__ == "__main__":
    unittest.main()
