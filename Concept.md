# Spatial Deconvolution Pipeline

Mapping single-cell RNA-seq reference data onto 10x Visium spatial spots to estimate
cell-type composition and spatial organization in tissue — built in two phases, going
from a validated analysis package to a containerized, orchestrated pipeline.

**Status:** 🚧 Phase 1 in progress

---

## Why two phases

Deconvolution modeling (getting DestVI/Tangram to converge and produce biologically
sensible results) and pipeline orchestration (Nextflow, Docker, reproducible execution)
are each non-trivial on their own. Building both at once risks a pipeline that runs but
whose outputs are unvalidated, or good science with no engineering around it. Phase 1
proves the analysis is correct in plain Python. Phase 2 wraps that proven logic in a
production-style workflow.

---

## Phase 1 — Analysis Package (Python)

**Goal:** a tested, documented Python package that performs the full analysis on a
public 10x Visium lymph node dataset + matched scRNA-seq reference, callable from
notebooks but with zero analysis logic living *in* the notebooks.

**Deliverables**
- [ ] `src/spatialdeconv/` package: `qc.py`, `reference.py`, `deconvolve.py`, `spatial_stats.py`, `plotting.py`
- [ ] `deconvolve.py` exposes one interface for both methods: `run_deconvolution(adata_sp, adata_ref, method="destvi"|"tangram")`
- [ ] Annotated scRNA-seq reference (cell types via scanpy/Leiden + marker genes)
- [ ] Per-spot cell-type proportions via DestVI (primary) and Tangram (baseline comparison)
- [ ] squidpy spatial stats: neighborhood enrichment, co-occurrence by dominant cell type
- [ ] `tests/` running the full pipeline on a subsampled fixture (~200 spots / ~500 cells) in under a minute
- [ ] GitHub Actions CI running the test suite on every push
- [ ] `environment.yml` with pinned versions
- [ ] Notebooks: `01_reference_annotation.ipynb`, `02_deconvolution.ipynb`, `03_spatial_patterns.ipynb`
- [ ] README hero figure: tissue map with cell-type overlay

**Tech stack:** Python, Scanpy, AnnData, scvi-tools (DestVI), Tangram, squidpy, pytest

**Definition of done:** a stranger clones the repo, runs 3 commands, and reproduces the
hero figure — plus can articulate the DestVI vs. Tangram comparison in one paragraph.

### Step-by-step build order

**Step 1 — Repo skeleton + environment**
```bash
mkdir spatial-deconv-project && cd spatial-deconv-project
git init
mkdir -p src/spatialdeconv data/raw data/processed notebooks tests/fixtures configs .github/workflows
touch environment.yml README.md
```
Write `environment.yml` (python=3.10, scanpy, squidpy, scvi-tools, anndata, pytorch, pytest,
pip: tangram-sc). Create the env, activate it, confirm `import scanpy, squidpy, scvi` works
before writing any analysis code.

**Step 2 — Get the data**
Write `data/download_data.py`: downloads the 10x Visium lymph node dataset and the matched
scRNA-seq reference into `data/raw/` (both are linked from 10x Genomics' public datasets page
and from the scvi-tools/squidpy tutorial pages — search "10x Visium human lymph node dataset"
and "scvi-tools DestVI tutorial data" to find the exact URLs). Add a checksum check. Run it once,
confirm both `.h5ad`/`.h5` files load in a scratch notebook.

**Step 3 — QC (`src/spatialdeconv/qc.py`)**
Write `qc_spatial(adata, min_counts, min_genes, pct_mito)` and `qc_reference(adata, ...)`.
Test each on the real data in a scratch notebook first (look at the distributions, pick
sensible thresholds), *then* move the finalized logic into the `.py` file. Don't write the
function blind — eyeball the QC plots first.

**Step 4 — Reference annotation (`src/spatialdeconv/reference.py`)**
Normalize → log-transform → HVG selection → PCA → neighbors → Leiden clustering on the
scRNA-seq reference. Look up marker genes for expected lymph node cell types (B cells,
T cells, macrophages, etc.) and label clusters by marker expression. Save the annotated
object with a `cell_type` column. Sanity check: does a UMAP colored by `cell_type` look
like clean, separated clusters?

**Step 5 — Deconvolution baseline first (`src/spatialdeconv/deconvolve.py`)**
Implement **Tangram first**, not DestVI — it's simpler and faster, so it proves your data
loading, alignment between spatial/reference gene sets, and output format all work before
you take on DestVI's extra complexity. Get `run_deconvolution(..., method="tangram")`
producing a proportions matrix you can plot.

**Step 6 — Add DestVI as the primary method**
Follow the scvi-tools DestVI tutorial structure: train `CondSCVI` on the reference, then
`DestVI` on the spatial data, extract proportions into the *same* output shape Tangram
produced in Step 5. This is where the shared interface pays off — same downstream plotting
code works for both. Expect this step to take the longest; model training + debugging
convergence is the real work here.

**Step 7 — Spatial stats (`src/spatialdeconv/spatial_stats.py`)**
Using the DestVI proportions, assign each spot a dominant cell type. Build a spatial
neighbors graph with `sq.gr.spatial_neighbors`, run `sq.gr.nhood_enrichment` and a
co-occurrence score. Plot which cell types spatially cluster together.

**Step 8 — Notebooks**
Now write the three notebooks — they should only *call* your package functions and produce
figures/narrative, no analysis logic. If a notebook cell has more than a few lines of
non-plotting code, that logic belongs in the package instead.

**Step 9 — Tests + CI**
Subsample the real data to ~200 spots / ~500 cells, save as `tests/fixtures/`. Write tests
that run each pipeline stage on this fixture and check output shapes/types (not biological
correctness — just "does it run and produce the right structure"). Add
`.github/workflows/ci.yml` running `pytest` on every push.

**Step 10 — README polish**
Add the hero figure (tissue map with cell-type overlay), a 3-command reproduce section, and
a short paragraph stating the biological question and what you found.

---

## Phase 2 — Orchestrated Pipeline (Nextflow + Docker)

**Goal:** refactor the validated Phase 1 steps into a Nextflow workflow, so the project
demonstrates pipeline engineering on top of already-proven analysis logic — modeled on
production pipelines like `scpca-nf`.

**Deliverables**
- [ ] `main.nf` + `modules/` — one Nextflow process per Phase 1 pipeline stage (QC → reference annotation → deconvolution → spatial stats)
- [ ] `Dockerfile` packaging all dependencies (or per-process containers if scope grows)
- [ ] `nextflow.config` with a `test` profile that runs on the same small fixture data used in Phase 1 CI
- [ ] Architecture diagram in the README: data flow from input matrices → deconvolution → spatial map
- [ ] CI: pipeline runs end-to-end on fixture data in a clean container on every push
- [ ] Decision documented explicitly: pipeline starts from **processed count matrices**, not raw FASTQ (alignment/counting is out of scope — noted as a possible Phase 3)

**Tech stack:** Nextflow, Docker, (Phase 1 stack, containerized)

**Definition of done:** `nextflow run main.nf -profile test,docker` reproduces Phase 1's
output from a clean environment with no manual setup.

### Step-by-step build order

**Step 1 — Install Nextflow, run a hello-world pipeline**
Before touching your real pipeline, install Nextflow and run the official "hello world"
example to confirm the tool itself works on your machine. Don't debug Nextflow syntax and
your own pipeline logic at the same time.

**Step 2 — One process at a time, starting from the simplest stage**
Convert Phase 1's QC step into a single Nextflow process first (input: raw `.h5ad`, output:
QC'd `.h5ad`, running your Phase 1 `qc.py` functions inside a small wrapper script). Get
`nextflow run` executing just this one process successfully before adding the next.

**Step 3 — Chain the remaining stages**
Add reference annotation, then deconvolution, then spatial stats as their own processes,
wiring outputs of one into inputs of the next via Nextflow channels. Test after each
addition — a pipeline that fails is much easier to debug with one new process at a time
than four at once.

**Step 4 — Dockerize**
Write a `Dockerfile` that installs your `environment.yml` dependencies. Build it locally,
confirm your package imports correctly inside the container, then point each Nextflow
process at that image via `nextflow.config`.

**Step 5 — Test profile**
Add a `test` profile in `nextflow.config` that points at the small fixture data from Phase 1
instead of the full dataset, so the whole pipeline runs in under a few minutes for CI and
for anyone trying the repo without downloading gigabytes of data.

**Step 6 — CI**
Add a GitHub Actions job that installs Nextflow + Docker and runs
`nextflow run main.nf -profile test,docker`, failing the build if it doesn't complete.

**Step 7 — Architecture diagram + README**
Draw the data flow (input matrices → QC → reference annotation → deconvolution → spatial
stats → output map) as a simple diagram (excalidraw, mermaid, or even a clean hand-drawn
export) and add it near the top of the README, above the Phase 1 hero figure.

---

## Data

- 10x Visium human/mouse lymph node (spatial)
- Matched scRNA-seq lymph node reference (single-cell)

Both public, no auth/paywall/ethics approval required. Fetched via `data/download_data.py`,
not committed to git.

## Repository structure (target, end of Phase 2)

```
spatial-deconv-project/
  data/                  download_data.py
  src/spatialdeconv/      qc.py, reference.py, deconvolve.py, spatial_stats.py, plotting.py
  notebooks/               01_reference_annotation.ipynb, 02_deconvolution.ipynb, 03_spatial_patterns.ipynb
  main.nf                  # Phase 2
  modules/                 # Phase 2: one .nf file per pipeline stage
  Dockerfile                # Phase 2
  nextflow.config           # Phase 2
  tests/fixtures/           tiny_visium.h5ad, tiny_ref.h5ad
  environment.yml
  .github/workflows/ci.yml
  README.md
```

## Roadmap

- [x] Scope defined
- [ ] Phase 1: analysis package + CI
- [ ] Phase 2: Nextflow + Docker
- [ ] Phase 3 (stretch): raw FASTQ → count matrix step; second tissue type (e.g. tumor Visium) for generalization
