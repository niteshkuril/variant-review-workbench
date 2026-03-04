# Deployment Notes

This document describes the current hosted shape for the Variant Review Workbench web app.

## Scope

The hosted service is a research-use demo over the existing pipeline. It is not a clinical system, a PHI-grade platform, or a multi-tenant privacy boundary.

## Runtime Shape

- web framework: Flask
- application entrypoint: `src.web.app`
- job execution: in-process threaded jobs by default
- persistent storage requirement: yes
- persistent cache target: ClinVar processed SQLite cache

## Render Baseline

Recommended Render service settings for the current repository state:

- Root Directory: blank
- Build Command: `pip install -r requirements.txt`
- Start Command: `flask --app src.web.app run --host 0.0.0.0 --port $PORT`
- Instance Type: `Starter`
- Disk Mount Path: `/var/data`

The start command still uses the Flask development server because `gunicorn` is not currently included in the project dependencies.

## Persistent Disk Layout

Recommended mounted layout:

- uploads: `/var/data/uploads`
- run outputs: `/var/data/runs`
- ClinVar raw files: `/var/data/clinvar/raw`
- ClinVar processed cache: `/var/data/clinvar/processed/clinvar_lookup_cache.sqlite3`

Raw ClinVar files expected by the app:

- `/var/data/clinvar/raw/variant_summary.txt.gz`
- `/var/data/clinvar/raw/summary_of_conflicting_interpretations.txt`
- `/var/data/clinvar/raw/submission_summary.txt.gz`

## Environment Variables

The web app reads deployment settings from these environment variables:

- `VRW_JOB_EXECUTION_MODE`
- `VRW_MAX_UPLOAD_MB`
- `VRW_UPLOAD_ROOT`
- `VRW_RUN_OUTPUT_ROOT`
- `VRW_RUN_RETENTION_HOURS`
- `VRW_CLINVAR_RAW_DIR`
- `VRW_CLINVAR_PROCESSED_DIR`
- `VRW_CLINVAR_VARIANT_SUMMARY`
- `VRW_CLINVAR_CONFLICT_SUMMARY`
- `VRW_CLINVAR_SUBMISSION_SUMMARY`
- `VRW_CLINVAR_CACHE_DB`
- `VRW_DISABLE_CLINVAR_CACHE`
- optional `VRW_DATA_ROOT`

Suggested Render values:

```text
VRW_JOB_EXECUTION_MODE=threaded
VRW_MAX_UPLOAD_MB=25
VRW_UPLOAD_ROOT=/var/data/uploads
VRW_RUN_OUTPUT_ROOT=/var/data/runs
VRW_RUN_RETENTION_HOURS=24
VRW_CLINVAR_RAW_DIR=/var/data/clinvar/raw
VRW_CLINVAR_PROCESSED_DIR=/var/data/clinvar/processed
VRW_CLINVAR_VARIANT_SUMMARY=/var/data/clinvar/raw/variant_summary.txt.gz
VRW_CLINVAR_CONFLICT_SUMMARY=/var/data/clinvar/raw/summary_of_conflicting_interpretations.txt
VRW_CLINVAR_SUBMISSION_SUMMARY=/var/data/clinvar/raw/submission_summary.txt.gz
VRW_CLINVAR_CACHE_DB=/var/data/clinvar/processed/clinvar_lookup_cache.sqlite3
VRW_DISABLE_CLINVAR_CACHE=false
```

## Health Check Behavior

The app exposes `/healthz`.

- `200 ok` means configured storage roots and required ClinVar paths are present
- `503 degraded` means one or more required reference paths are missing

For first deployment on a fresh disk, `/healthz` will stay degraded until the ClinVar files exist at the configured mounted paths. If Render needs a temporary health check during initial bootstrap, `/` can be used until the raw files are uploaded.

## First Hosted Run

The first hosted run against a fresh raw ClinVar snapshot is expected to be slow.

What happens:

- the app builds the processed ClinVar SQLite cache on the persistent disk
- that preprocessing cost is paid once per raw snapshot state
- later runs reuse the warmed cache and should be much faster

In practical terms, hosted users do not each pay a cold-start preprocessing penalty once the cache already exists on the persistent disk.

## Local Web Startup

Install dependencies:

```powershell
python -m pip install -r requirements.txt
```

Run the web app locally:

```powershell
flask --app src.web.app run --host 127.0.0.1 --port 5000
```

If you want the local web app to use an isolated runtime workspace instead of repository-relative defaults, set `VRW_DATA_ROOT` before launching.

## Privacy And Retention

- uploads are stored in per-run directories
- stale workspaces are removed opportunistically on new submissions based on `VRW_RUN_RETENTION_HOURS`
- uploaded files can persist until another submission triggers cleanup
- do not upload protected health information to the hosted demo

## Smoke Test Checklist

After deployment, verify:

1. `GET /` renders the upload form.
2. `GET /docs` renders the project overview.
3. `GET /healthz` returns `200` once ClinVar files are present.
4. A VCF submission reaches the results page.
5. `report`, `html`, `json`, and `md` outputs are reachable for a successful run.
