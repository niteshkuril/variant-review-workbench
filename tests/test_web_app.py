from __future__ import annotations

import csv
import gzip
import io
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from src.web import create_app


class WebAppTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = tempfile.TemporaryDirectory()
        root = Path(self.tmpdir.name)
        clinvar_root = root / "clinvar"
        self.variant_summary = clinvar_root / "variant_summary.txt.gz"
        self.conflict_summary = clinvar_root / "summary_of_conflicting_interpretations.txt"
        self.submission_summary = clinvar_root / "submission_summary.txt.gz"
        self.cache_db = root / "clinvar_lookup_cache.sqlite3"
        clinvar_root.mkdir(parents=True, exist_ok=True)
        self._write_reference_files()
        self.app = create_app(
            {
                "TESTING": True,
                "JOB_EXECUTION_MODE": "inline",
                "UPLOAD_ROOT": str(root / "uploads"),
                "RUN_OUTPUT_ROOT": str(root / "runs"),
                "RUN_RETENTION_HOURS": 1,
                "CLINVAR_VARIANT_SUMMARY": str(self.variant_summary),
                "CLINVAR_CONFLICT_SUMMARY": str(self.conflict_summary),
                "CLINVAR_SUBMISSION_SUMMARY": str(self.submission_summary),
                "CLINVAR_CACHE_DB": str(self.cache_db),
                "DISABLE_CLINVAR_CACHE": False,
            }
        )
        self.client = self.app.test_client()

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def _write_reference_files(self) -> None:
        with gzip.open(self.variant_summary, "wt", encoding="utf-8", newline="") as handle:
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
                    "TP53 example",
                    "TP53",
                    "Pathogenic",
                    "Jan 01, 2025",
                    "Li-Fraumeni syndrome",
                    "reviewed by expert panel",
                    "germline",
                    "GRCh38",
                    "17",
                    "1234",
                    "43045702",
                    "A",
                    "G",
                    "RCV000000001",
                ]
            )

        with self.conflict_summary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                [
                    "#Gene_Symbol",
                    "NCBI_Variation_ID",
                    "ClinVar_Preferred",
                    "Submitter1",
                    "Submitter1_SCV",
                    "Submitter1_ClinSig",
                    "Submitter1_LastEval",
                    "Submitter1_ReviewStatus",
                    "Submitter1_Sub_Condition",
                    "Submitter1_Description",
                    "Submitter2",
                    "Submitter2_SCV",
                    "Submitter2_ClinSig",
                    "Submitter2_LastEval",
                    "Submitter2_ReviewStatus",
                    "Submitter2_Sub_Condition",
                    "Submitter2_Description",
                    "Rank_diff",
                    "Conflict_Reported",
                    "Variant_type",
                    "Submitter1_Method",
                    "Submitter2_Method",
                ]
            )
            writer.writerow(
                [
                    "TP53",
                    "1234",
                    "TP53 example",
                    "Lab A",
                    "SCV1",
                    "Pathogenic",
                    "Jan 01, 2025",
                    "criteria provided, single submitter",
                    "",
                    "",
                    "Lab B",
                    "SCV2",
                    "Uncertain significance",
                    "Jan 02, 2025",
                    "criteria provided, single submitter",
                    "",
                    "",
                    "1",
                    "yes",
                    "SNV",
                    "clinical testing",
                    "clinical testing",
                ]
            )

        with gzip.open(self.submission_summary, "wt", encoding="utf-8", newline="") as handle:
            handle.write("##Comment\n")
            handle.write(
                "#VariationID\tClinicalSignificance\tDateLastEvaluated\tDescription\t"
                "SubmittedPhenotypeInfo\tReportedPhenotypeInfo\tReviewStatus\tCollectionMethod\t"
                "OriginCounts\tSubmitter\tSCV\tSubmittedGeneSymbol\tExplanationOfInterpretation\t"
                "SomaticClinicalImpact\tOncogenicity\tContributesToAggregateClassification\n"
            )
            handle.write(
                "1234\tPathogenic\tJan 01, 2025\t-\tLi-Fraumeni syndrome\t"
                "C0085390:Li-Fraumeni syndrome\treviewed by expert panel\tclinical testing\t"
                "germline:1\tLab A\tSCV1\tTP53\t-\t-\t-\tyes\n"
            )

    @staticmethod
    def _demo_vcf_bytes() -> bytes:
        return (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "17\t43045702\t.\tA\tG\t100\tPASS\tGENE=TP53;IMPACT=HIGH\n"
        ).encode("utf-8")

    @classmethod
    def _demo_vcfgz_bytes(cls) -> bytes:
        buffer = io.BytesIO()
        with gzip.GzipFile(fileobj=buffer, mode="wb") as handle:
            handle.write(cls._demo_vcf_bytes())
        return buffer.getvalue()

    def test_home_page_renders_form_shell(self) -> None:
        response = self.client.get("/")

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Review a VCF through the same pipeline used by the CLI.", response.data)
        self.assertIn(b"Run Setup", response.data)
        self.assertIn(b'enctype="multipart/form-data"', response.data)
        self.assertIn(b"3 Formats", response.data)
        self.assertIn(b"Research-use only.", response.data)
        self.assertIn(b"25 MB", response.data)
        self.assertIn(b"Individual uploads do not trigger a full cache build.", response.data)

    def test_docs_page_renders_project_context(self) -> None:
        response = self.client.get("/docs")

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Variant Review Workbench", response.data)
        self.assertIn(b"Browser workflow", response.data)
        self.assertIn(b"Hosted demo guardrails", response.data)
        self.assertIn(b"processed ClinVar cache should be built offline", response.data)

    def test_health_check_returns_ok(self) -> None:
        response = self.client.get("/healthz")

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        assert payload is not None
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["job_execution_mode"], "inline")
        self.assertEqual(payload["max_upload_mb"], 25)
        self.assertTrue(payload["paths"]["upload_root"].endswith("uploads"))
        self.assertTrue(payload["paths"]["run_output_root"].endswith("runs"))
        self.assertTrue(payload["checks"]["clinvar_variant_summary_exists"])

    def test_health_check_returns_degraded_for_missing_reference_paths(self) -> None:
        broken_app = create_app(
            {
                "TESTING": True,
                "JOB_EXECUTION_MODE": "inline",
                "UPLOAD_ROOT": str(Path(self.tmpdir.name) / "broken_uploads"),
                "RUN_OUTPUT_ROOT": str(Path(self.tmpdir.name) / "broken_runs"),
                "CLINVAR_VARIANT_SUMMARY": str(Path(self.tmpdir.name) / "missing" / "variant_summary.txt.gz"),
                "CLINVAR_CONFLICT_SUMMARY": str(Path(self.tmpdir.name) / "missing" / "summary_of_conflicting_interpretations.txt"),
                "CLINVAR_SUBMISSION_SUMMARY": str(Path(self.tmpdir.name) / "missing" / "submission_summary.txt.gz"),
                "CLINVAR_CACHE_DB": str(Path(self.tmpdir.name) / "missing_processed" / "clinvar_lookup_cache.sqlite3"),
            }
        )

        response = broken_app.test_client().get("/healthz")

        self.assertEqual(response.status_code, 503)
        payload = response.get_json()
        assert payload is not None
        self.assertEqual(payload["status"], "degraded")
        self.assertFalse(payload["checks"]["clinvar_variant_summary_exists"])
        self.assertFalse(payload["checks"]["clinvar_conflict_summary_exists"])
        self.assertFalse(payload["checks"]["clinvar_submission_summary_exists"])

    def test_create_run_redirects_to_results_shell(self) -> None:
        response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh38",
                "export_format": "json",
                "vcf_file": (io.BytesIO(self._demo_vcf_bytes()), "demo.vcf"),
            },
            content_type="multipart/form-data",
        )

        self.assertEqual(response.status_code, 302)
        self.assertIn("/runs/run-", response.headers["Location"])

    def test_export_only_submission_redirects_with_export_preference(self) -> None:
        response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh38",
                "mode": "export_only",
                "export_format": "md",
                "vcf_file": (io.BytesIO(self._demo_vcfgz_bytes()), "demo.vcf.gz"),
            },
            content_type="multipart/form-data",
        )

        self.assertEqual(response.status_code, 302)
        self.assertIn("/runs/export-", response.headers["Location"])
        self.assertIn("/export/md", response.headers["Location"])

        redirected = self.client.get(response.headers["Location"])
        self.assertEqual(redirected.status_code, 200)
        self.assertIn("text/markdown", redirected.content_type)
        self.assertIn(b"# Variant Review Report", redirected.data)

    def test_results_page_renders_placeholder(self) -> None:
        create_response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh38",
                "export_format": "html",
                "vcf_file": (io.BytesIO(self._demo_vcf_bytes()), "demo.vcf"),
            },
            content_type="multipart/form-data",
        )
        response = self.client.get(create_response.headers["Location"])

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Review the generated report and exported artifacts.", response.data)
        self.assertIn(b"status-badge succeeded", response.data)
        self.assertIn(b"html", response.data)
        self.assertIn(b"The report and export artifacts are ready.", response.data)
        self.assertIn(b"Uploaded VCF", response.data)
        self.assertIn(b"Report Preview", response.data)
        self.assertIn(b"Export HTML", response.data)
        self.assertIn(b"Export JSON", response.data)
        self.assertIn(b"Export Markdown", response.data)
        self.assertIn(b"Non-clinical output only.", response.data)

    def test_results_page_shows_no_match_warning_when_clinvar_match_count_is_zero(self) -> None:
        create_response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh37",
                "export_format": "html",
                "vcf_file": (io.BytesIO(self._demo_vcf_bytes()), "demo.vcf"),
            },
            content_type="multipart/form-data",
        )
        response = self.client.get(create_response.headers["Location"])

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"No ClinVar matches were found for this run.", response.data)
        self.assertIn(b"All variants were assigned", response.data)
        self.assertIn(b"context_only", response.data)
        self.assertIn(b"GRCh37", response.data)

    def test_status_endpoint_returns_job_result(self) -> None:
        create_response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh37",
                "export_format": "json",
                "vcf_file": (io.BytesIO(self._demo_vcf_bytes()), "demo.vcf"),
            },
            content_type="multipart/form-data",
        )
        run_path = create_response.headers["Location"]
        run_id = run_path.rstrip("/").split("/")[-1]

        response = self.client.get(f"/runs/{run_id}/status")

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        assert payload is not None
        self.assertEqual(payload["job_id"], run_id)
        self.assertEqual(payload["status"], "succeeded")
        self.assertEqual(payload["metadata"]["assembly"], "GRCh37")
        self.assertEqual(payload["result"]["summary"]["clinvar_matched_count"], 0)
        self.assertEqual(payload["result"]["summary"]["clinvar_unmatched_count"], 1)

    def test_status_endpoint_returns_pipeline_result(self) -> None:
        create_response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh38",
                "export_format": "json",
                "vcf_file": (io.BytesIO(self._demo_vcf_bytes()), "demo.vcf"),
            },
            content_type="multipart/form-data",
        )
        run_id = create_response.headers["Location"].rstrip("/").split("/")[-1]

        response = self.client.get(f"/runs/{run_id}/status")

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        assert payload is not None
        self.assertEqual(payload["job_id"], run_id)
        self.assertEqual(payload["status"], "succeeded")
        self.assertEqual(payload["metadata"]["assembly"], "GRCh38")
        self.assertEqual(payload["result"]["mode"], "report")
        self.assertTrue(payload["metadata"]["uploaded_vcf_path"].endswith("demo.vcf"))
        self.assertIn(run_id, payload["metadata"]["output_dir"])
        self.assertTrue(payload["result"]["report_html_path"].endswith("report.html"))
        self.assertEqual(payload["result"]["summary"]["clinvar_matched_count"], 1)

    def test_status_endpoint_returns_404_for_missing_run(self) -> None:
        response = self.client.get("/runs/run-missing/status")

        self.assertEqual(response.status_code, 404)
        self.assertEqual(response.get_json(), {"error": "run not found"})

    def test_report_route_returns_404_for_missing_run(self) -> None:
        response = self.client.get("/runs/run-missing/report")

        self.assertEqual(response.status_code, 404)

    def test_export_route_returns_404_for_missing_run(self) -> None:
        response = self.client.get("/runs/run-missing/export/json")

        self.assertEqual(response.status_code, 404)

    def test_create_run_rejects_missing_upload(self) -> None:
        response = self.client.post(
            "/runs",
            data={"assembly": "GRCh38", "export_format": "json"},
            content_type="multipart/form-data",
        )

        self.assertEqual(response.status_code, 400)
        self.assertIn(b"Submission error", response.data)
        self.assertIn(b"A VCF or VCF.GZ file is required.", response.data)

    def test_create_run_rejects_non_vcf_upload(self) -> None:
        response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh38",
                "export_format": "json",
                "vcf_file": (io.BytesIO(b"not a vcf"), "notes.txt"),
            },
            content_type="multipart/form-data",
        )

        self.assertEqual(response.status_code, 400)
        self.assertIn(b"Uploaded file must end with .vcf or .vcf.gz.", response.data)

    def test_create_run_rejects_non_vcf_contents_with_vcf_extension(self) -> None:
        response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh38",
                "export_format": "json",
                "vcf_file": (io.BytesIO(b"plain text"), "notes.vcf"),
            },
            content_type="multipart/form-data",
        )

        self.assertEqual(response.status_code, 400)
        self.assertIn(b"Uploaded file does not appear to be a VCF.", response.data)

    def test_create_run_rejects_invalid_gzip_payload(self) -> None:
        response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh38",
                "export_format": "json",
                "vcf_file": (io.BytesIO(b"not gzip data"), "notes.vcf.gz"),
            },
            content_type="multipart/form-data",
        )

        self.assertEqual(response.status_code, 400)
        self.assertIn(b"Uploaded .vcf.gz file is not gzip-compressed.", response.data)

    def test_create_run_rejects_request_over_upload_limit(self) -> None:
        small_limit_app = create_app(
            {
                "TESTING": True,
                "JOB_EXECUTION_MODE": "inline",
                "MAX_UPLOAD_MB": 1,
                "MAX_CONTENT_LENGTH": 1024 * 1024,
                "UPLOAD_ROOT": str(Path(self.tmpdir.name) / "small_uploads"),
                "RUN_OUTPUT_ROOT": str(Path(self.tmpdir.name) / "small_runs"),
                "CLINVAR_VARIANT_SUMMARY": str(self.variant_summary),
                "CLINVAR_CONFLICT_SUMMARY": str(self.conflict_summary),
                "CLINVAR_SUBMISSION_SUMMARY": str(self.submission_summary),
                "CLINVAR_CACHE_DB": str(self.cache_db),
                "DISABLE_CLINVAR_CACHE": False,
            }
        )

        response = small_limit_app.test_client().post(
            "/runs",
            data={
                "assembly": "GRCh38",
                "export_format": "json",
                "vcf_file": (io.BytesIO(self._demo_vcf_bytes() + (b"A" * (1024 * 1024))), "demo.vcf"),
            },
            content_type="multipart/form-data",
        )

        self.assertEqual(response.status_code, 413)
        self.assertIn(b"Uploaded file exceeds the 1 MB limit.", response.data)

    def test_uploaded_file_is_isolated_in_run_workspace(self) -> None:
        response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh38",
                "export_format": "json",
                "vcf_file": (io.BytesIO(self._demo_vcf_bytes()), "..\\unsafe-demo.vcf"),
            },
            content_type="multipart/form-data",
        )
        run_id = response.headers["Location"].rstrip("/").split("/")[-1]

        status_response = self.client.get(f"/runs/{run_id}/status")
        payload = status_response.get_json()
        assert payload is not None

        uploaded_path = Path(payload["metadata"]["uploaded_vcf_path"])
        output_dir = Path(payload["metadata"]["output_dir"])
        self.assertTrue(uploaded_path.exists())
        self.assertEqual(uploaded_path.parent.name, run_id)
        self.assertEqual(output_dir.parent.name, run_id)
        self.assertEqual(uploaded_path.name, "unsafe-demo.vcf")

    def test_report_route_serves_generated_html(self) -> None:
        create_response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh38",
                "export_format": "html",
                "vcf_file": (io.BytesIO(self._demo_vcf_bytes()), "demo.vcf"),
            },
            content_type="multipart/form-data",
        )
        run_id = create_response.headers["Location"].rstrip("/").split("/")[-1]

        response = self.client.get(f"/runs/{run_id}/report")

        self.assertEqual(response.status_code, 200)
        self.assertIn("text/html", response.content_type)
        self.assertIn(b"Variant Review Report", response.data)

    def test_export_route_serves_generated_json_and_markdown(self) -> None:
        create_response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh38",
                "export_format": "json",
                "vcf_file": (io.BytesIO(self._demo_vcf_bytes()), "demo.vcf"),
            },
            content_type="multipart/form-data",
        )
        run_id = create_response.headers["Location"].rstrip("/").split("/")[-1]

        json_response = self.client.get(f"/runs/{run_id}/export/json")
        markdown_response = self.client.get(f"/runs/{run_id}/export/md")

        self.assertEqual(json_response.status_code, 200)
        self.assertIn("application/json", json_response.content_type)
        self.assertIn(b'"report_title": "Variant Review Report"', json_response.data)
        self.assertIn(b'"summary"', json_response.data)

        self.assertEqual(markdown_response.status_code, 200)
        self.assertIn("text/markdown", markdown_response.content_type)
        self.assertIn(b"# Variant Review Report", markdown_response.data)
        self.assertIn(b"## Top Findings", markdown_response.data)

    def test_export_route_returns_404_for_unsupported_format(self) -> None:
        create_response = self.client.post(
            "/runs",
            data={
                "assembly": "GRCh38",
                "export_format": "json",
                "vcf_file": (io.BytesIO(self._demo_vcf_bytes()), "demo.vcf"),
            },
            content_type="multipart/form-data",
        )
        run_id = create_response.headers["Location"].rstrip("/").split("/")[-1]

        response = self.client.get(f"/runs/{run_id}/export/xml")

        self.assertEqual(response.status_code, 404)

    def test_threaded_failed_job_surfaces_status_and_results_page(self) -> None:
        threaded_app = create_app(
            {
                "TESTING": True,
                "JOB_EXECUTION_MODE": "threaded",
                "UPLOAD_ROOT": str(Path(self.tmpdir.name) / "threaded_uploads"),
                "RUN_OUTPUT_ROOT": str(Path(self.tmpdir.name) / "threaded_runs"),
                "RUN_RETENTION_HOURS": 1,
                "CLINVAR_VARIANT_SUMMARY": str(self.variant_summary),
                "CLINVAR_CONFLICT_SUMMARY": str(self.conflict_summary),
                "CLINVAR_SUBMISSION_SUMMARY": str(self.submission_summary),
                "CLINVAR_CACHE_DB": str(self.cache_db),
                "DISABLE_CLINVAR_CACHE": False,
            }
        )
        client = threaded_app.test_client()

        with patch("src.web.app.run_pipeline_with_result", side_effect=RuntimeError("synthetic pipeline failure")):
            create_response = client.post(
                "/runs",
                data={
                    "assembly": "GRCh38",
                    "export_format": "json",
                    "vcf_file": (io.BytesIO(self._demo_vcf_bytes()), "demo.vcf"),
                },
                content_type="multipart/form-data",
            )

        self.assertEqual(create_response.status_code, 302)
        run_id = create_response.headers["Location"].rstrip("/").split("/")[-1]

        status_payload = None
        for _ in range(100):
            status_response = client.get(f"/runs/{run_id}/status")
            status_payload = status_response.get_json()
            assert status_payload is not None
            if status_payload["status"] == "failed":
                break
            time.sleep(0.02)

        assert status_payload is not None
        self.assertEqual(status_payload["status"], "failed")
        self.assertEqual(status_payload["error"], "synthetic pipeline failure")

        results_response = client.get(f"/runs/{run_id}")
        self.assertEqual(results_response.status_code, 200)
        self.assertIn(b"Run failed", results_response.data)
        self.assertIn(b"synthetic pipeline failure", results_response.data)


if __name__ == "__main__":
    unittest.main()
