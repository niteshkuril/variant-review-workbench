"""In-process job execution primitives for the web interface."""

from __future__ import annotations

from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import UTC, datetime
import json
from pathlib import Path
from threading import Lock
from typing import Any, Callable


def _utcnow_iso() -> str:
    """Return a stable UTC timestamp string for job metadata."""
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


@dataclass(slots=True)
class JobRecord:
    """Represents one submitted web job."""

    job_id: str
    status: str
    created_at: str
    updated_at: str
    mode: str
    export_format: str | None
    result: dict[str, Any] | None = None
    error: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


class JobStore:
    """Thread-safe in-memory store for submitted web jobs."""

    def __init__(self, *, state_root: Path | None = None) -> None:
        self._jobs: dict[str, JobRecord] = {}
        self._lock = Lock()
        self._state_root = state_root

    @staticmethod
    def _record_to_payload(job: JobRecord) -> dict[str, Any]:
        return {
            "job_id": job.job_id,
            "status": job.status,
            "created_at": job.created_at,
            "updated_at": job.updated_at,
            "mode": job.mode,
            "export_format": job.export_format,
            "result": job.result,
            "error": job.error,
            "metadata": job.metadata,
        }

    @staticmethod
    def _record_from_payload(payload: dict[str, Any]) -> JobRecord:
        return JobRecord(
            job_id=str(payload.get("job_id", "")),
            status=str(payload.get("status", "queued")),
            created_at=str(payload.get("created_at", _utcnow_iso())),
            updated_at=str(payload.get("updated_at", _utcnow_iso())),
            mode=str(payload.get("mode", "report")),
            export_format=payload.get("export_format"),
            result=payload.get("result"),
            error=payload.get("error"),
            metadata=dict(payload.get("metadata") or {}),
        )

    def _state_path(self, job_id: str) -> Path | None:
        if self._state_root is None:
            return None
        return self._state_root / job_id / "job_state.json"

    def _persist_job_state(self, job: JobRecord) -> None:
        state_path = self._state_path(job.job_id)
        if state_path is None:
            return

        state_path.parent.mkdir(parents=True, exist_ok=True)
        payload = self._record_to_payload(job)
        temp_path = state_path.with_suffix(state_path.suffix + ".tmp")
        temp_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        temp_path.replace(state_path)

    def _load_job_state(self, job_id: str) -> JobRecord | None:
        state_path = self._state_path(job_id)
        if state_path is None or not state_path.exists():
            return None

        try:
            payload = json.loads(state_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        if not isinstance(payload, dict):
            return None
        return self._record_from_payload(payload)

    def create_job(self, *, job_id: str, mode: str, export_format: str | None, metadata: dict[str, Any]) -> JobRecord:
        """Create and store a new queued job record."""
        timestamp = _utcnow_iso()
        job = JobRecord(
            job_id=job_id,
            status="queued",
            created_at=timestamp,
            updated_at=timestamp,
            mode=mode,
            export_format=export_format,
            metadata=dict(metadata),
        )
        with self._lock:
            self._jobs[job_id] = job
            self._persist_job_state(job)
        return job

    def get_job(self, job_id: str) -> JobRecord | None:
        """Return a job record by identifier."""
        with self._lock:
            existing = self._jobs.get(job_id)
            if existing is not None:
                return existing
            recovered = self._load_job_state(job_id)
            if recovered is not None:
                self._jobs[job_id] = recovered
            return recovered

    def start_job(self, job_id: str) -> None:
        """Mark a job as running."""
        with self._lock:
            job = self._jobs[job_id]
            job.status = "running"
            job.updated_at = _utcnow_iso()
            self._persist_job_state(job)

    def complete_job(self, job_id: str, result: dict[str, Any]) -> None:
        """Mark a job as successfully completed."""
        with self._lock:
            job = self._jobs[job_id]
            job.status = "succeeded"
            job.result = result
            job.updated_at = _utcnow_iso()
            self._persist_job_state(job)

    def fail_job(self, job_id: str, error: str) -> None:
        """Mark a job as failed."""
        with self._lock:
            job = self._jobs[job_id]
            job.status = "failed"
            job.error = error
            job.updated_at = _utcnow_iso()
            self._persist_job_state(job)


class JobRunner:
    """Submit and track background jobs for the web interface."""

    def __init__(self, *, store: JobStore, execution_mode: str = "threaded", max_workers: int = 2) -> None:
        if execution_mode not in {"threaded", "inline"}:
            raise ValueError("execution_mode must be 'threaded' or 'inline'")

        self.store = store
        self.execution_mode = execution_mode
        self._executor = ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="vrw-web") if execution_mode == "threaded" else None

    def submit(
        self,
        *,
        job_id: str,
        mode: str,
        export_format: str | None,
        metadata: dict[str, Any],
        task: Callable[[], dict[str, Any]],
    ) -> JobRecord:
        """Submit a new job and return its initial record."""
        job = self.store.create_job(job_id=job_id, mode=mode, export_format=export_format, metadata=metadata)

        if self.execution_mode == "inline":
            self._run_job(job_id, task)
            return job

        assert self._executor is not None
        future = self._executor.submit(self._run_job, job_id, task)
        future.add_done_callback(self._consume_future_exception)
        return job

    def _run_job(self, job_id: str, task: Callable[[], dict[str, Any]]) -> None:
        self.store.start_job(job_id)
        try:
            result = task()
        except Exception as error:  # pragma: no cover - exercised through public API
            self.store.fail_job(job_id, str(error))
            raise
        self.store.complete_job(job_id, result)

    @staticmethod
    def _consume_future_exception(future: Future[None]) -> None:
        """Drain worker exceptions to avoid unobserved-future warnings."""
        try:
            future.result()
        except Exception:
            return
