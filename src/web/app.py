"""Web application entrypoint for variant-review-workbench."""

from __future__ import annotations

import argparse
from pathlib import Path
from uuid import uuid4

from flask import Flask, Response, abort, jsonify, redirect, render_template, request, send_file, url_for
from werkzeug.exceptions import RequestEntityTooLarge

from ..app_service import PipelineUsageError, run_pipeline_with_result
from ..models import GenomeAssembly
from ..report_builder import write_markdown_report, write_report_export_json
from .jobs import JobRunner, JobStore
from .settings import WebRuntimeSettings
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
    pipeline_result = run_pipeline_with_result(pipeline_args)
    markdown_path = write_markdown_report(output_dir / "report.md", pipeline_result.report_context)
    report_json_path = write_report_export_json(output_dir / "report.json", pipeline_result.report_context)
    return {
        "job_id": job_id,
        "assembly": assembly,
        "mode": mode,
        "export_format": export_format,
        "pharmgkb_enabled": pharmgkb_enabled,
        "uploaded_vcf_path": str(uploaded_vcf_path),
        "output_dir": str(output_dir),
        "report_html_path": str(pipeline_result.outputs["report_html"]),
        "report_markdown_path": str(markdown_path),
        "report_json_path": str(report_json_path),
        "summary_json_path": str(pipeline_result.outputs["summary_json"]),
        "prioritized_variants_json_path": str(pipeline_result.outputs["prioritized_variants_json"]),
        "annotated_variants_csv_path": str(pipeline_result.outputs["annotated_variants_csv"]),
        "run_metadata_json_path": str(pipeline_result.outputs["run_metadata_json"]),
        "summary": pipeline_result.run_metadata.statistics.model_dump(mode="json"),
        "message": "Pipeline completed successfully and the generated report is available below.",
    }


def _resolve_export_path(job_result: dict[str, object], export_format: str) -> tuple[Path, str, str]:
    """Resolve a supported export format to a file path, mimetype, and download name."""
    mapping = {
        "html": ("report_html_path", "text/html; charset=utf-8", "report.html"),
        "md": ("report_markdown_path", "text/markdown; charset=utf-8", "report.md"),
        "json": ("report_json_path", "application/json; charset=utf-8", "report.json"),
    }
    if export_format not in mapping:
        raise PipelineUsageError("export format must be one of: html, md, json")

    field_name, mimetype, download_name = mapping[export_format]
    export_path = Path(str(job_result[field_name]))
    return export_path, mimetype, download_name


def create_app(test_config: dict | None = None) -> Flask:
    """Create the Flask application for the thin web interface."""
    project_root = Path(__file__).resolve().parents[2]
    template_dir = project_root / "templates"
    static_dir = Path(__file__).resolve().with_name("static")

    app = Flask(__name__, template_folder=str(template_dir), static_folder=str(static_dir))
    runtime_settings = WebRuntimeSettings.from_env(project_root)
    app.config.from_mapping(runtime_settings.to_flask_config())
    if test_config:
        app.config.update(test_config)

    runtime_settings = WebRuntimeSettings(
        project_root=project_root,
        job_execution_mode=str(app.config["JOB_EXECUTION_MODE"]),
        max_upload_mb=int(app.config["MAX_UPLOAD_MB"]),
        job_max_workers=max(1, int(app.config["JOB_MAX_WORKERS"])),
        upload_root=Path(str(app.config["UPLOAD_ROOT"])),
        run_output_root=Path(str(app.config["RUN_OUTPUT_ROOT"])),
        run_retention_hours=int(app.config["RUN_RETENTION_HOURS"]),
        clinvar_variant_summary=Path(str(app.config["CLINVAR_VARIANT_SUMMARY"])),
        clinvar_conflict_summary=Path(str(app.config["CLINVAR_CONFLICT_SUMMARY"])) if app.config["CLINVAR_CONFLICT_SUMMARY"] else None,
        clinvar_submission_summary=Path(str(app.config["CLINVAR_SUBMISSION_SUMMARY"])) if app.config["CLINVAR_SUBMISSION_SUMMARY"] else None,
        clinvar_cache_db=Path(str(app.config["CLINVAR_CACHE_DB"])) if app.config["CLINVAR_CACHE_DB"] else None,
        disable_clinvar_cache=bool(app.config["DISABLE_CLINVAR_CACHE"]),
    )
    app.extensions["runtime_settings"] = runtime_settings

    app.config["JOB_EXECUTION_MODE_CONFIGURED"] = str(app.config["JOB_EXECUTION_MODE"])

    upload_root = runtime_settings.upload_root
    run_output_root = runtime_settings.run_output_root
    ensure_storage_roots(upload_root, run_output_root)

    app.extensions["job_store"] = JobStore(state_root=run_output_root)
    app.extensions["job_runner"] = JobRunner(
        store=app.extensions["job_store"],
        execution_mode=str(app.config["JOB_EXECUTION_MODE"]),
        max_workers=max(1, int(app.config["JOB_MAX_WORKERS"])),
    )

    @app.context_processor
    def inject_web_scope() -> dict[str, object]:
        return {
            "web_scope": {
                "max_upload_mb": int(app.config["MAX_UPLOAD_MB"]),
                "run_retention_hours": int(app.config["RUN_RETENTION_HOURS"]),
            }
        }

    @app.errorhandler(RequestEntityTooLarge)
    def handle_request_too_large(_: RequestEntityTooLarge) -> tuple[str, int]:
        return (
            render_template(
                "web/home.html.j2",
                page_title="Variant Review Workbench",
                current_page="home",
                form_error=f"Uploaded file exceeds the {app.config['MAX_UPLOAD_MB']} MB limit.",
            ),
            413,
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
        if mode == "export_only" and app.config["JOB_EXECUTION_MODE"] == "inline":
            return redirect(url_for("run_export", run_id=job_id, export_format=export_format))
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
        return send_file(report_path, mimetype="text/html; charset=utf-8")

    @app.get("/runs/<run_id>/export/<export_format>")
    def run_export(run_id: str, export_format: str) -> Response:
        job = app.extensions["job_store"].get_job(run_id)
        if job is None or job.status != "succeeded" or not job.result:
            abort(404)

        try:
            export_path, mimetype, download_name = _resolve_export_path(job.result, export_format)
        except PipelineUsageError:
            abort(404)

        if not export_path.exists():
            abort(404)

        response = Response(export_path.read_text(encoding="utf-8"), mimetype=mimetype)
        response.headers["Content-Disposition"] = f'attachment; filename="{download_name}"'
        return response

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
    def healthz() -> tuple[dict[str, object], int]:
        payload = runtime_settings.health_snapshot()
        payload["job_execution_mode"] = str(app.config["JOB_EXECUTION_MODE"])
        payload["job_execution_mode_configured"] = str(app.config["JOB_EXECUTION_MODE_CONFIGURED"])
        payload["job_max_workers"] = int(app.config["JOB_MAX_WORKERS"])
        payload["max_upload_mb"] = runtime_settings.max_upload_mb
        status_code = 200 if payload["status"] == "ok" else 503
        return payload, status_code

    return app


app = create_app()
