"""Runtime settings for the web interface and hosted deployments."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


def _env_flag(name: str, default: bool) -> bool:
    """Parse a boolean environment variable."""
    raw = os.getenv(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def _env_int(name: str, default: int) -> int:
    """Parse an integer environment variable."""
    raw = os.getenv(name)
    if raw is None:
        return default
    return int(raw)


@dataclass(slots=True)
class WebRuntimeSettings:
    """Resolved storage and reference-data settings for the web app."""

    project_root: Path
    job_execution_mode: str
    max_upload_mb: int
    upload_root: Path
    run_output_root: Path
    run_retention_hours: int
    clinvar_variant_summary: Path
    clinvar_conflict_summary: Path | None
    clinvar_submission_summary: Path | None
    clinvar_cache_db: Path | None
    disable_clinvar_cache: bool

    @classmethod
    def from_env(cls, project_root: Path) -> "WebRuntimeSettings":
        """Resolve settings from environment variables with local defaults."""
        data_root = Path(os.getenv("VRW_DATA_ROOT", str(project_root / ".web_runtime")))
        clinvar_raw_dir = Path(os.getenv("VRW_CLINVAR_RAW_DIR", str(project_root / "data" / "clinvar" / "raw")))
        clinvar_processed_dir = Path(os.getenv("VRW_CLINVAR_PROCESSED_DIR", str(project_root / "data" / "clinvar" / "processed")))

        return cls(
            project_root=project_root,
            job_execution_mode=os.getenv("VRW_JOB_EXECUTION_MODE", "threaded"),
            max_upload_mb=_env_int("VRW_MAX_UPLOAD_MB", 25),
            upload_root=Path(os.getenv("VRW_UPLOAD_ROOT", str(data_root / "uploads"))),
            run_output_root=Path(os.getenv("VRW_RUN_OUTPUT_ROOT", str(data_root / "runs"))),
            run_retention_hours=_env_int("VRW_RUN_RETENTION_HOURS", 24),
            clinvar_variant_summary=Path(
                os.getenv("VRW_CLINVAR_VARIANT_SUMMARY", str(clinvar_raw_dir / "variant_summary.txt.gz"))
            ),
            clinvar_conflict_summary=Path(
                os.getenv(
                    "VRW_CLINVAR_CONFLICT_SUMMARY",
                    str(clinvar_raw_dir / "summary_of_conflicting_interpretations.txt"),
                )
            ),
            clinvar_submission_summary=Path(
                os.getenv("VRW_CLINVAR_SUBMISSION_SUMMARY", str(clinvar_raw_dir / "submission_summary.txt.gz"))
            ),
            clinvar_cache_db=Path(
                os.getenv("VRW_CLINVAR_CACHE_DB", str(clinvar_processed_dir / "clinvar_lookup_cache.sqlite3"))
            ),
            disable_clinvar_cache=_env_flag("VRW_DISABLE_CLINVAR_CACHE", False),
        )

    def to_flask_config(self) -> dict[str, object]:
        """Serialize settings into the Flask config shape used by the app."""
        return {
            "JOB_EXECUTION_MODE": self.job_execution_mode,
            "MAX_CONTENT_LENGTH": self.max_upload_mb * 1024 * 1024,
            "MAX_UPLOAD_MB": self.max_upload_mb,
            "UPLOAD_ROOT": str(self.upload_root),
            "RUN_OUTPUT_ROOT": str(self.run_output_root),
            "RUN_RETENTION_HOURS": self.run_retention_hours,
            "CLINVAR_VARIANT_SUMMARY": str(self.clinvar_variant_summary),
            "CLINVAR_CONFLICT_SUMMARY": str(self.clinvar_conflict_summary) if self.clinvar_conflict_summary else None,
            "CLINVAR_SUBMISSION_SUMMARY": str(self.clinvar_submission_summary) if self.clinvar_submission_summary else None,
            "CLINVAR_CACHE_DB": str(self.clinvar_cache_db) if self.clinvar_cache_db else None,
            "DISABLE_CLINVAR_CACHE": self.disable_clinvar_cache,
        }

    def health_snapshot(self) -> dict[str, object]:
        """Build a deployment-oriented health snapshot for configured paths."""
        cache_parent_exists = True
        if self.clinvar_cache_db is not None:
            cache_parent_exists = self.clinvar_cache_db.parent.exists()

        checks = {
            "upload_root_exists": self.upload_root.exists(),
            "run_output_root_exists": self.run_output_root.exists(),
            "clinvar_variant_summary_exists": self.clinvar_variant_summary.exists(),
            "clinvar_conflict_summary_exists": (
                self.clinvar_conflict_summary.exists() if self.clinvar_conflict_summary is not None else True
            ),
            "clinvar_submission_summary_exists": (
                self.clinvar_submission_summary.exists() if self.clinvar_submission_summary is not None else True
            ),
            "clinvar_cache_parent_exists": cache_parent_exists,
        }
        status = "ok" if all(checks.values()) else "degraded"
        return {
            "status": status,
            "checks": checks,
            "paths": {
                "upload_root": str(self.upload_root),
                "run_output_root": str(self.run_output_root),
                "clinvar_variant_summary": str(self.clinvar_variant_summary),
                "clinvar_conflict_summary": str(self.clinvar_conflict_summary) if self.clinvar_conflict_summary else None,
                "clinvar_submission_summary": str(self.clinvar_submission_summary) if self.clinvar_submission_summary else None,
                "clinvar_cache_db": str(self.clinvar_cache_db) if self.clinvar_cache_db else None,
            },
        }
