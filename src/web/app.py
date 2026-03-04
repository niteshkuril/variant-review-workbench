"""Web application entrypoint for variant-review-workbench."""

from __future__ import annotations

import argparse
from pathlib import Path
from uuid import uuid4

from flask import Flask, Response, abort, jsonify, redirect, render_template, request, url_for

from ..app_service import PipelineUsageError, run_pipeline_with_details
from ..models import GenomeAssembly
from .jobs import JobRunner, JobStore
from .storage import (
    UploadValidationError,
    cleanup_expired_run_directories,
    create_run_workspace,
    ensure_storage_roots,
    save_uploaded_vcf,
)


def _parse_assembly(value: str) -> GenomeAssembly:
    """Parse a web assembly form value into a supported enum."""
    normalized = value.strip()
    if normalized == GenomeAssembly.GRCH37.value:
        return GenomeAssembly.GRCH37
    if normalized == GenomeAssembly.GRCH38.value:
        return GenomeAssembly.GRCH38
    raise PipelineUsageError("assembly must be either 'GRCh37' or 'GRCh38'")


def _build_pipeline_task_args(
    *,
    uploaded_vcf_path: Path,
    output_dir: Path,
    assembly: str,
    pharmgkb_enabled: bool,
    app: Flask,
) -> argparse.Namespace:
    """Build the shared pipeline argument namespace for a web-submitted run."""
    return argparse.Namespace(
        input=str(uploaded_vcf_path),
        assembly=_parse_assembly(assembly),
        variant_summary=app.config["CLINVAR_VARIANT_SUMMARY"],
        conflict_summary=app.config["CLINVAR_CONFLICT_SUMMARY"],
        submission_summary=app.config["CLINVAR_SUBMISSION_SUMMARY"],
        clinvar_cache_db=app.config["CLINVAR_CACHE_DB"],
        disable_clinvar_cache=bool(app.config["DISABLE_CLINVAR_CACHE"]),
        out_dir=str(output_dir),
        enable_pharmgkb=pharmgkb_enabled,
    )


def _run_pipeline_job(
    *,
    job_id: str,
    assembly: str,
    mode: str,
    pharmgkb_enabled: bool,
    export_format: str | None,
    uploaded_vcf_path: Path,
    output_dir: Path,
    app: Flask,
) -> dict[str, object]:
    """Execute the shared reporting pipeline for one web-submitted run."""
    pipeline_args = _build_pipeline_task_args(
        uploaded_vcf_path=uploaded_vcf_path,
        output_dir=output_dir,
        assembly=assembly,
        pharmgkb_enabled=pharmgkb_enabled,
        app=app,
    )
    outputs, run_metadata = run_pipeline_with_details(pipeline_args)
    return {
        "job_id": job_id,
        "assembly": assembly,
        "mode": mode,
        "export_format": export_format,
        "pharmgkb_enabled": pharmgkb_enabled,
        "uploaded_vcf_path": str(uploaded_vcf_path),
        "output_dir": str(output_dir),
        "report_html_path": str(outputs["report_html"]),
        "summary_json_path": str(outputs["summary_json"]),
        "prioritized_variants_json_path": str(outputs["prioritized_variants_json"]),
        "annotated_variants_csv_path": str(outputs["annotated_variants_csv"]),
        "run_metadata_json_path": str(outputs["run_metadata_json"]),
        "summary": run_metadata.statistics.model_dump(mode="json"),
        "message": "Pipeline completed successfully and the generated report is available below.",
    }


def create_app(test_config: dict | None = None) -> Flask:
    """Create the Flask application for the thin web interface."""
    project_root = Path(__file__).resolve().parents[2]
    template_dir = project_root / "templates"
    static_dir = Path(__file__).resolve().with_name("static")

    app = Flask(__name__, template_folder=str(template_dir), static_folder=str(static_dir))
    app.config["MAX_CONTENT_LENGTH"] = 25 * 1024 * 1024
    app.config["JOB_EXECUTION_MODE"] = "threaded"
    app.config["UPLOAD_ROOT"] = str(project_root / ".web_runtime" / "uploads")
    app.config["RUN_OUTPUT_ROOT"] = str(project_root / ".web_runtime" / "runs")
    app.config["RUN_RETENTION_HOURS"] = 24
    app.config["CLINVAR_VARIANT_SUMMARY"] = str(project_root / "data" / "clinvar" / "raw" / "variant_summary.txt.gz")
    app.config["CLINVAR_CONFLICT_SUMMARY"] = str(project_root / "data" / "clinvar" / "raw" / "summary_of_conflicting_interpretations.txt")
    app.config["CLINVAR_SUBMISSION_SUMMARY"] = str(project_root / "data" / "clinvar" / "raw" / "submission_summary.txt.gz")
    app.config["CLINVAR_CACHE_DB"] = str(project_root / "data" / "clinvar" / "processed" / "clinvar_lookup_cache.sqlite3")
    app.config["DISABLE_CLINVAR_CACHE"] = False
    if test_config:
        app.config.update(test_config)

    upload_root = Path(app.config["UPLOAD_ROOT"])
    run_output_root = Path(app.config["RUN_OUTPUT_ROOT"])
    ensure_storage_roots(upload_root, run_output_root)

    app.extensions["job_store"] = JobStore()
    app.extensions["job_runner"] = JobRunner(
        store=app.extensions["job_store"],
        execution_mode=app.config["JOB_EXECUTION_MODE"],
    )

    @app.get("/")
    def home() -> str:
        return render_template(
            "web/home.html.j2",
            page_title="Variant Review Workbench",
            current_page="home",
        )

    @app.get("/docs")
    def docs() -> str:
        return render_template(
            "web/docs.html.j2",
            page_title="Workbench Docs",
            current_page="docs",
        )

    @app.post("/runs")
    def create_run() -> tuple[str, int] | object:
        mode = "export_only" if request.form.get("mode") == "export_only" else "report"
        job_prefix = "export" if mode == "export_only" else "run"
        job_id = f"{job_prefix}-{uuid4().hex[:8]}"
        export_format = request.form.get("export_format", "json")
        assembly = request.form.get("assembly", "GRCh38")
        pharmgkb_enabled = request.form.get("enable_pharmgkb") == "true"
        uploaded_vcf = request.files.get("vcf_file")

        cleanup_expired_run_directories(
            upload_root=upload_root,
            run_output_root=run_output_root,
            retention_hours=int(app.config["RUN_RETENTION_HOURS"]),
        )
        try:
            workspace = create_run_workspace(job_id=job_id, upload_root=upload_root, run_output_root=run_output_root)
            workspace = save_uploaded_vcf(upload=uploaded_vcf, workspace=workspace)
        except UploadValidationError as error:
            return (
                render_template(
                    "web/home.html.j2",
                    page_title="Variant Review Workbench",
                    current_page="home",
                    form_error=str(error),
                ),
                400,
            )

        metadata = {
            "assembly": assembly,
            "pharmgkb_enabled": pharmgkb_enabled,
            "requested_export_format": export_format,
            "uploaded_vcf_path": str(workspace.uploaded_vcf_path),
            "output_dir": str(workspace.output_dir),
            "cleanup_policy": f"Run workspaces older than {app.config['RUN_RETENTION_HOURS']} hour(s) are removed opportunistically on new submissions.",
        }

        job_runner: JobRunner = app.extensions["job_runner"]
        job_runner.submit(
            job_id=job_id,
            mode=mode,
            export_format=export_format,
            metadata=metadata,
            task=lambda: _run_pipeline_job(
                job_id=job_id,
                assembly=assembly,
                mode=mode,
                pharmgkb_enabled=pharmgkb_enabled,
                export_format=export_format,
                uploaded_vcf_path=workspace.uploaded_vcf_path,
                output_dir=workspace.output_dir,
                app=app,
            ),
        )
        return redirect(url_for("results", run_id=job_id))

    @app.get("/runs/<run_id>")
    def results(run_id: str) -> str:
        job = app.extensions["job_store"].get_job(run_id)
        if job is None:
            abort(404)

        return render_template(
            "web/results.html.j2",
            page_title="Run Results",
            current_page="results",
            run_id=run_id,
            export_format=job.export_format,
            job=job,
        )

    @app.get("/runs/<run_id>/report")
    def run_report(run_id: str) -> Response:
        job = app.extensions["job_store"].get_job(run_id)
        if job is None or job.status != "succeeded" or not job.result:
            abort(404)

        report_path = Path(str(job.result["report_html_path"]))
        if not report_path.exists():
            abort(404)
        return Response(report_path.read_text(encoding="utf-8"), mimetype="text/html")

    @app.get("/runs/<run_id>/status")
    def run_status(run_id: str) -> tuple[object, int]:
        job = app.extensions["job_store"].get_job(run_id)
        if job is None:
            return {"error": "run not found"}, 404

        return jsonify(
            {
                "job_id": job.job_id,
                "status": job.status,
                "mode": job.mode,
                "export_format": job.export_format,
                "created_at": job.created_at,
                "updated_at": job.updated_at,
                "metadata": job.metadata,
                "result": job.result,
                "error": job.error,
            }
        ), 200

    @app.get("/healthz")
    def healthz() -> tuple[dict[str, str], int]:
        return {
            "status": "ok",
            "job_execution_mode": app.config["JOB_EXECUTION_MODE"],
            "upload_root": str(upload_root),
            "run_output_root": str(run_output_root),
        }, 200

    return app


app = create_app()
