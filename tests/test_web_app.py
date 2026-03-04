from __future__ import annotations

import csv
import gzip
import io
import tempfile
import unittest
from pathlib import Path

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
        self.assertIn(b"Run the workbench without touching the command line.", response.data)
        self.assertIn(b"Run Setup", response.data)
        self.assertIn(b'enctype="multipart/form-data"', response.data)

    def test_docs_page_renders_project_context(self) -> None:
        response = self.client.get("/docs")

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Variant Review Workbench", response.data)
        self.assertIn(b"thin analyst-facing layer", response.data)

    def test_health_check_returns_ok(self) -> None:
        response = self.client.get("/healthz")

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        assert payload is not None
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["job_execution_mode"], "inline")
        self.assertTrue(payload["upload_root"].endswith("uploads"))
        self.assertTrue(payload["run_output_root"].endswith("runs"))

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
        self.assertIn(b"Results Shell", response.data)
        self.assertIn(b"Status:", response.data)
        self.assertIn(b"html", response.data)
        self.assertIn(b"Pipeline completed successfully", response.data)
        self.assertIn(b"Uploaded VCF:", response.data)
        self.assertIn(b"Report Preview", response.data)
        self.assertIn(b"Export HTML", response.data)
        self.assertIn(b"Export JSON", response.data)
        self.assertIn(b"Export Markdown", response.data)

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


if __name__ == "__main__":
    unittest.main()
