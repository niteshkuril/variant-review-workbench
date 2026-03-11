# Variant Review Workbench

ClinVar-first small-variant triage and reporting tool for research-oriented review workflows.

The workbench accepts a VCF, matches variants against a local ClinVar snapshot, highlights conflicting interpretations, optionally enriches findings with PharmGKB, ranks the review queue, and emits analyst-friendly HTML and machine-readable outputs.

## Problem Statement

VCFs are compact and machine-friendly, but they are not ideal review artifacts. Analysts often need to answer a narrower question first:

- Which variants matched a known ClinVar record?
- Which findings have conflicting interpretations?
- Which records should be reviewed first?
- Which variants may have pharmacogenomics context worth surfacing?

This repository focuses on that gap. It is not a full annotation platform or a clinical interpretation engine. It is a reproducible, inspectable triage workbench for small-variant review.

## What The Tool Does

- reads `.vcf` and `.vcf.gz` inputs
- normalizes one record per alternate allele
- matches variants to a local ClinVar snapshot using assembly-aware coordinate and allele keys
- attaches conflict and submission context when available
- optionally enriches variants with PharmGKB gene, variant, clinical annotation, and guideline data
- ranks findings with transparent heuristics
- writes HTML, JSON, CSV, and run metadata outputs

## Product Boundary

This project is:

- a research triage workbench
- a reproducible annotation and reporting tool

This project is not:

- a clinical decision support system
- an ACMG classifier
- a star-allele caller
- a treatment recommendation engine

## Architecture

```text
VCF / VCF.GZ
    |
    v
VCF parser
    |
    v
Normalized InputVariant records
    |
    v
ClinVar exact-match index
    |
    +--> conflict attachment
    |
    +--> submission evidence attachment
    |
    v
AnnotatedVariant records
    |
    +--> optional PharmGKB enrichment
    |
    v
RankedVariant records
    |
    +--> HTML report
    +--> prioritized_variants.json
    +--> annotated_variants.csv
    +--> summary.json
    `--> run_metadata.json
```

## Repository Layout

```text
variant-review-workbench/
|-- src/
|   |-- annotator.py
|   |-- cli.py
|   |-- clinvar_index.py
|   |-- models.py
|   |-- pgx_enrichment.py
|   |-- ranker.py
|   |-- report_builder.py
|   `-- vcf_parser.py
|-- templates/
|   `-- report.html.j2
|-- data/
|   |-- clinvar/
|   |-- pharmgkb/
|   |-- references/
|   `-- demo.vcf
|-- tests/
|-- README.md
`-- pyproject.toml
```

## Inputs

### Required

- input VCF or VCF.GZ
- ClinVar `variant_summary.txt.gz`
- reference assembly: `GRCh37` or `GRCh38`

### Optional

- `summary_of_conflicting_interpretations.txt`
- `submission_summary.txt.gz`
- PharmGKB enrichment via `--enable-pharmgkb`

### Download Locations

Download the required ClinVar files from the official NCBI FTP `tab_delimited` directory:

- directory listing: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/`
- `variant_summary.txt.gz`: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz`
- `summary_of_conflicting_interpretations.txt`: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/summary_of_conflicting_interpretations.txt`
- `submission_summary.txt.gz`: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz`

Recommended local placement for the examples in this README:

- `data/clinvar/raw/variant_summary.txt.gz`
- `data/clinvar/raw/summary_of_conflicting_interpretations.txt`
- `data/clinvar/raw/submission_summary.txt.gz`

## Outputs

Each run writes:

- `annotated_variants.csv`
  - stable CSV export of ranked variant records
  - list-valued fields are serialized as JSON arrays inside each cell
- `prioritized_variants.json`
  - machine-readable prioritized variant artifact with `schema_version`, `artifact_type`, and `records`
- `summary.json`
  - machine-readable summary artifact with stable count fields and priority-tier counts
- `run_metadata.json`
  - reproducibility metadata, source provenance, and counts
- `report.html`
  - analyst-facing HTML report with top findings, conflicts, methods, and limitations

## Important Runtime Note

This tool now uses a persistent processed ClinVar cache.

- by default the CLI builds and reuses a SQLite cache at `data/clinvar/processed/clinvar_lookup_cache.sqlite3`
- the first local run against a new raw ClinVar snapshot is a preprocessing run and can take a long time
- after that cache exists, repeated runs against the same snapshot should be much faster

Observed timing on the staged repository data:

- first run with cache build: about `1402` seconds, about `23 minutes 22 seconds`
- warm run reusing the cache: about `2.83` seconds

What this means in practice:

- if a fresh local environment appears slow on the first real run, that is expected
- that first local run is building a reusable queryable index from the raw ClinVar files
- the warm-run path is the intended day-to-day workflow
- hosted deployments should not rely on the first public request to perform this build

Cache controls:

- default cache location: `data/clinvar/processed/clinvar_lookup_cache.sqlite3`
- override location: `--clinvar-cache-db <path>`
- disable cache and force raw-file reads: `--disable-clinvar-cache`
- prebuild or refresh the cache without running a report: `python -m src.cache_bootstrap ...`

## Output Contract

The machine-readable outputs are intentionally separate from the human-oriented HTML report.

- `prioritized_variants.json` is the canonical structured export for downstream scripts
- `annotated_variants.csv` uses the same field set as the JSON records where practical
- list-valued fields such as `condition_names`, `flags`, and `ranking_rationale` remain lists in JSON and are encoded as JSON arrays in CSV cells
- `summary.json` uses stable count names:
  - `input_variant_count`
  - `clinvar_matched_count`
  - `clinvar_unmatched_count`
  - `conflict_flagged_count`
  - `pharmgkb_enriched_count`
  - `review_priority_tier_counts`

## Setup

```powershell
python -m pip install -r requirements.txt
```

## Quick Start

Pick one of these entry points:

1. CLI: run the existing pipeline directly and inspect the generated artifacts in `outputs/`.
2. Web: start the Flask app locally and submit a VCF from the browser.

To inspect the generated outputs first, start with the web interface or open `outputs/demo_run/report.html`.

## Web Interface

The repository also includes a thin Flask web layer that reuses the same backend pipeline as the CLI.

Local startup:

```powershell
flask --app src.web.app run --host 127.0.0.1 --port 5000
```

Then open `http://127.0.0.1:5000`.

Web flow:

1. Upload a `.vcf` or `.vcf.gz`.
2. Select `GRCh37` or `GRCh38`.
3. Optionally enable PharmGKB enrichment.
4. Either open the browser report or go straight to `html`, `json`, or `md` export.

For deployment-specific details, see [DEPLOYMENT.md](DEPLOYMENT.md).

Hosted process summary:

- the web app reuses a shared processed ClinVar cache rather than building one per uploaded VCF
- `VRW_RUN_RETENTION_HOURS` applies to temporary run workspaces, not to cache warm-up or user cooldown
- the recommended hosted workflow is to build `clinvar_lookup_cache.sqlite3` offline and upload it to the mounted disk before public use

Hosted/runtime environment variables:

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

Recommended hosted disk layout:

- uploads: `/var/data/uploads`
- run outputs: `/var/data/runs`
- ClinVar raw files: `/var/data/clinvar/raw`
- ClinVar cache: `/var/data/clinvar/processed/clinvar_lookup_cache.sqlite3`

Hosted recommendation:

- do not rely on the first public web submission to build the shared ClinVar cache
- prebuild the processed cache offline, then place the finished SQLite file on the mounted disk before treating the site as ready

### Web Guardrails

- the web app is a convenience interface over the same research pipeline, not a clinical product
- upload acceptance is limited to `.vcf` and `.vcf.gz` files and the request size cap defaults to `25` MB
- the hosted demo stores submitted files in per-run workspaces and removes stale workspaces opportunistically after the configured retention window
- the hosted demo is not a PHI-grade privacy boundary, so protected health information should not be uploaded
- `/healthz` is intended to return healthy only after the configured ClinVar source files and cache parent directories exist on disk

## Example Usage

### Fastest Local Review Paths

CLI reviewer path:

```powershell
python -m src.cli `
  --input data\demo.vcf `
  --assembly GRCh38 `
  --variant-summary data\clinvar\raw\variant_summary.txt.gz `
  --conflict-summary data\clinvar\raw\summary_of_conflicting_interpretations.txt `
  --submission-summary data\clinvar\raw\submission_summary.txt.gz `
  --out-dir outputs\demo_run
```

Web reviewer path:

```powershell
flask --app src.web.app run --host 127.0.0.1 --port 5000
```

Then submit `data\demo.vcf` through the homepage.

### Base ClinVar Run

This command will build the processed ClinVar cache on first use if it does not already exist.

```powershell
python -m src.cli `
  --input data\demo.vcf `
  --assembly GRCh38 `
  --variant-summary data\clinvar\raw\variant_summary.txt.gz `
  --conflict-summary data\clinvar\raw\summary_of_conflicting_interpretations.txt `
  --submission-summary data\clinvar\raw\submission_summary.txt.gz `
  --out-dir outputs\demo_run
```

### Cache Control Examples

Use a specific cache path:

```powershell
python -m src.cli `
  --input data\demo.vcf `
  --assembly GRCh38 `
  --variant-summary data\clinvar\raw\variant_summary.txt.gz `
  --conflict-summary data\clinvar\raw\summary_of_conflicting_interpretations.txt `
  --submission-summary data\clinvar\raw\submission_summary.txt.gz `
  --clinvar-cache-db data\clinvar\processed\clinvar_lookup_cache.sqlite3 `
  --out-dir outputs\demo_run
```

Disable the processed cache entirely:

```powershell
python -m src.cli `
  --input data\demo.vcf `
  --assembly GRCh38 `
  --variant-summary data\clinvar\raw\variant_summary.txt.gz `
  --conflict-summary data\clinvar\raw\summary_of_conflicting_interpretations.txt `
  --submission-summary data\clinvar\raw\submission_summary.txt.gz `
  --disable-clinvar-cache `
  --out-dir outputs\demo_run_no_cache
```

Prebuild or refresh the processed ClinVar cache without running a report:

```powershell
python -m src.cache_bootstrap `
  --variant-summary data\clinvar\raw\variant_summary.txt.gz `
  --conflict-summary data\clinvar\raw\summary_of_conflicting_interpretations.txt `
  --submission-summary data\clinvar\raw\submission_summary.txt.gz `
  --clinvar-cache-db data\clinvar\processed\clinvar_lookup_cache.sqlite3 `
  --force-rebuild
```

### Run With PharmGKB Enrichment

```powershell
python -m src.cli `
  --input data\demo.vcf `
  --assembly GRCh38 `
  --variant-summary data\clinvar\raw\variant_summary.txt.gz `
  --conflict-summary data\clinvar\raw\summary_of_conflicting_interpretations.txt `
  --submission-summary data\clinvar\raw\submission_summary.txt.gz `
  --out-dir outputs\demo_run_pgx `
  --enable-pharmgkb
```

## Demo Dataset

The repository demo input is intentionally small but not arbitrary. It is designed to show different review situations in one run:

- `BRCA1`
  - strong pathogenic ClinVar match with conflict surfaced
- `APC`
  - strong pathogenic ClinVar match without conflict
- `DPYD`
  - drug-response-oriented ClinVar match that becomes especially useful in the PharmGKB-enabled run
- `TP53`
  - conflict-heavy variant that stays important but scores below the strongest pathogenic findings

The demo file lives at:

- `data/demo.vcf`

On the current staged ClinVar snapshot, the base demo run produces:

- `4` input variants
- `4` ClinVar matches
- `3` conflict-flagged findings
- `3` `high_review_priority` findings
- `1` `review` finding

The PharmGKB-enabled demo run adds cached public PGx context and, on the current demo set, enriches all `4` variants.

## Demo Walkthrough

Inspect these generated artifacts in order:

1. `outputs/demo_run/report.html`
2. `outputs/demo_run/summary.json`
3. `outputs/demo_run/prioritized_variants.json`
4. `outputs/demo_run_pgx/report.html`

What the base demo should immediately show:

- the hero summary confirms all four demo variants matched ClinVar
- the top findings section shows both conflict-flagged and non-conflict findings
- the conflict review queue is not empty
- the rationale text explains why `BRCA1`, `APC`, and `DPYD` outrank `TP53`

What the PharmGKB demo should immediately show:

- the same ClinVar-first queue remains intact
- optional PGx context adds signal without replacing the core ClinVar interpretation layer
- the output provenance records both ClinVar and PharmGKB sources

## Demo Screenshots

Base report hero:

![Base report hero](data/screenshots/report_hero.png)

This view gives the fastest high-level read on the run:

- `4` input variants
- `4` ClinVar matches
- `3` conflict-flagged findings
- `3` high-priority findings

Base report top findings:

![Base report top findings](data/screenshots/top_findings.png)

This view shows the ranked findings near the top of the report:

- `BRCA1` as a high-priority conflict-flagged finding
- `APC` as a strong non-conflict pathogenic finding
- `DPYD` as a drug-response-oriented finding that becomes more interesting in the PharmGKB-enabled run

Base report conflict review queue:

![Base report conflict review queue](data/screenshots/conflict_review_queue.png)

This view shows the conflict review queue section:

- the conflict review queue is populated
- the highest-friction findings are isolated into a reviewable section

Base report variant table:

![Base report variant table](data/screenshots/variant_table.png)

This view shows the denser analyst-facing table:

- tier, score, ClinVar significance, review status, conflict flag, and workflow flags appear together
- the row-level output mirrors the machine-readable exports

PharmGKB-enabled source provenance:

![PharmGKB-enabled source provenance](data/screenshots/source_provenance.png)

This view shows the provenance section after a PharmGKB-enabled run:

- ClinVar sources remain explicit
- PharmGKB API and local PharmGKB cache are both recorded
- the optional enrichment layer is visible without obscuring the ClinVar-first workflow

## Ranking Approach

Ranking is heuristic and intentionally transparent.

The current score uses:

- ClinVar clinical significance
- ClinVar review strength
- conflict surfacing
- input impact severity
- optional PharmGKB context
- gene-symbol mismatch penalty

The system emits a numeric score, a priority tier, and a rationale list for every ranked variant.

Priority tiers:

- `high_review_priority`
- `review`
- `context_only`

## Example Review Questions This Tool Helps Answer

- Which variants have strong ClinVar support and should be reviewed first?
- Which findings are conflict-flagged and require closer inspection?
- Which unmatched records remain context-only?
- Which variants have optional PGx context worth surfacing for downstream review?

## Data Sources

### ClinVar

Used as the core local reference layer for variant matching and conflict attachment.

- FTP root: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/`
- maintenance and release notes: `https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/`
- FTP primer: `https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/`

Primary files used by this workbench:

- [`variant_summary.txt.gz`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)
- [`summary_of_conflicting_interpretations.txt`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/summary_of_conflicting_interpretations.txt)
- [`submission_summary.txt.gz`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz)

### PharmGKB / ClinPGx

Used only as optional enrichment.

- API base: `https://api.pharmgkb.org/v1`
- docs: `https://api.pharmgkb.org/`

Current integration uses stable public queries for:

- gene lookup by symbol
- variant lookup by symbol
- clinical annotations by gene symbol
- guideline annotations by gene symbol

### Source Reference Index

Human-readable source links are also maintained in:

- [external_sources.md](data/references/external_sources.md)

## Licensing And Attribution

This repository contains code written for the workbench itself. Upstream datasets and APIs remain governed by their respective providers.

- Repository code license: see [LICENSE](LICENSE)
  - free use for personal, educational, research, and internal non-commercial purposes
  - no commercial use
  - redistribution requires clear credit to the original author and project
- ClinVar data usage and redistribution expectations should be reviewed through NCBI documentation.
- PharmGKB API usage should follow PharmGKB terms and public API guidance.

Users of this repository should verify current upstream licensing and attribution requirements before redistributing derived datasets or packaging source snapshots.

## Testing

Run the current unit suite:

```powershell
python -m unittest discover -s tests -v
```

GitHub Actions also runs syntax checks and the unit suite on push and pull request.

The implemented system currently has coverage for:

- VCF parsing
- ClinVar index loading
- conflict and submission attachment
- annotation behavior
- ranking behavior
- report generation
- CLI orchestration
- PharmGKB caching, failure handling, and integration
- web route behavior
- web job execution and failure handling

## Current Status

Implemented and tested:

- ClinVar-first local annotation pipeline
- HTML report generation
- CSV and JSON exports
- optional PharmGKB enrichment
- end-to-end CLI orchestration
- thin Flask web interface with homepage, docs page, per-run uploads, report embedding, and exports

Current automated test count:

- `76` passing unit tests


## Limitations

- matching is exact and assembly-aware, but does not yet perform deeper variant normalization beyond the current key strategy
- unmatched variants are intentionally left as context-only rather than force-interpreted
- PharmGKB enrichment is optional and network-dependent when enabled
- this is a focused small-variant review tool, not a full-scale annotation framework

## Disclaimer

This tool is for research triage and educational review only. It is not intended for diagnosis, treatment selection, or other clinical decision-making.

The hosted web demo is also not a production privacy boundary. Do not submit protected health information or rely on uploaded-run retention as a formal data-governance control.
