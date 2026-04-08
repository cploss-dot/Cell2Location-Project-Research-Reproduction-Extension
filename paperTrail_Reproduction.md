# Cell2location Reproduction — Team 7 Paper Trail
**DS 596 Special Topics: Learning from Large-scale Biological Data**
**Team:** Vaidehi Gupta, Chloe Ploss, Addison Yam
**Professor:** Pawel Przytycki
**Paper:** Kleshchevnikov et al. (2022) "Cell2location maps fine-grained cell types in spatial transcriptomics" — Nature Biotechnology
**Goal:** Reproduce Figure 2b–f from the mouse brain analysis

---

## Environment Setup

### Cluster
- HPC cluster: SCC2 (Boston University)
- Working directory: `/projectnb/ds596/projects/Team 7/Cell2Location-Project-Research-Reproduction-Extention`
- Job scheduler: `qsub` (SGE/PBS-style). Note: `sbatch` (SLURM) is NOT available on SCC2 — BU uses `qsub` instead
- Singularity ONLY works inside a desktop interactive session, not a regular terminal login

### Only Working Workflow to Launch Jupyter
```bash
# Step 1: Open a desktop interactive session via the web portal
# Step 2: In the terminal inside that desktop session:
singularity shell cell2location.sif
# Step 3: Inside the singularity shell:
jupyter notebook
```
Everything must run in one continuous Jupyter notebook. No separate `.py` scripts or shell
scripts — variables like `adata_snrna` only exist in notebook memory and are inaccessible
to external scripts.

### Data Downloaded
**snRNA-seq (E-MTAB-11115):**
```bash
wget -r --no-parent ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/115/E-MTAB-11115/Files/
```
6 samples, ~2.8G total. Files used: `*_filtered_feature_bc_matrix.h5` and `cell_annotation.csv`.

**Visium spatial (E-MTAB-11114):**
```bash
wget -r --no-parent ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/114/E-MTAB-11114/Files/
```
5 sections, ~5.5G total. Files used: `*_filtered_feature_bc_matrix.h5` and extracted spatial folders.

### Visium Spatial Folder Extraction
```bash
cd mouseVisium/
for f in *.tar.gz; do
    sample="${f%%_spatial.tar.gz}"
    mkdir -p "${sample}_spatial"
    tar -xzf "$f" -C "${sample}_spatial"
done
```
Each extracted folder contains: `aligned_fiducials.jpg`, `detected_tissue_image.jpg`,
`scalefactors_json.json`, `tissue_hires_image.png`, `tissue_lowres_image.png`,
`tissue_positions_list.csv`.

**Important path note:** extraction creates `ST8059048_spatial/spatial/` — one extra
nesting level that affected how `sc.read_visium` was called (see difficulty below).

### Singularity Container
```bash
singularity pull cell2location.sif docker://quay.io/vitkl/cell2location:latest
```
Successfully downloaded ~4.8GB. Moved to data root directory.

### Container Version
The pulled container uses an **older cell2location API** — function-based rather than
class-based. This affected all downstream code (see difficulties below).

Confirmed working imports inside container:
- `cell2location` (no `__version__` attribute in this build — harmless)
- `scanpy`
- `torch`
- `bbknn`
- GPU: NOT available in interactive desktop session — all training runs on CPU

---

## Data Facts Confirmed
| Item | Value | Matches Paper? |
|------|-------|---------------|
| snRNA-seq samples | 6 | ✓ |
| Total cells before filtering | 44,429 | ✓ |
| Annotated cells after filtering | 40,531 | ✓ |
| Unique cell types (annotation_1) | 59 | ✓ |
| Visium sections | 5 (ST8059048–ST8059052) | ✓ |
| Cell ID format | `{sampleID}_{barcode}-1` | ✓ |

---

## Notebook Progress

### Cell 1 — Imports and paths ✅
```python
import scanpy as sc
import pandas as pd
import numpy as np
import os

data_dir = "data/"
scrna_dir = os.path.join(data_dir, "mousescRNAseq")
visium_dir = os.path.join(data_dir, "mouseVisium")
```

### Cell 2 — Load and concatenate all 6 snRNA-seq samples ✅
Loaded each sample with `sc.read_10x_h5`, prefixed barcodes with sample ID to avoid
collisions, concatenated with `sc.concat`.
- Output: `(44429, 31053)` combined

### Cell 3 — Load annotation file ✅
- `cell_annotation.csv`: 40,532 rows (including header), 4 columns
- 59 unique cell types confirmed via `annotation_1` column

### Cell 4 — Confirm cell ID format match ✅
Both loaded cell IDs and annotation IDs use `{sampleID}_{barcode}-1` format — confirmed
matching before merge.

### Cell 5 — Merge annotations and filter to annotated cells ✅
- Saved raw counts to `adata_snrna.layers["counts"]` before any normalization
- Merged `annotation_1` and `annotation_1_print` via index intersection
- Filtered to annotated cells only
- Output: ~40,531 cells, 59 cell types

### Cell 6 — Broad class mapping, normalize, PCA, BBKNN, UMAP ✅
- Defined 7 broad classes: Astrocytes, Excitatory, Inhibitory, Oligo-OPC, Micro, Endo, Other
- Normalized for visualization only (raw counts preserved in `layers["counts"]`)
- HVG selection: top 3,000 genes with `batch_key="sample"`
- PCA: 50 components
- Batch correction: BBKNN with `neighbors_within_batch=3` (matches paper)
- UMAP computed

### Cell 7 — Plot Figure 2b ✅
- Left panel: 59 fine-grained subtypes colored by `annotation_1`
- Right panel: broad class colored by `broad_class`
- Saved to `figure_2b.png`

### Cell 8 — Load all 5 Visium sections ✅
```python
import cell2location

visium_samples = ["ST8059048", "ST8059049", "ST8059050", "ST8059051", "ST8059052"]

slides = {}
for sample in visium_samples:
    h5_path = os.path.join(visium_dir, f"{sample}_filtered_feature_bc_matrix.h5")
    spatial_path = os.path.join(visium_dir, f"{sample}_spatial")
    adata_vis = sc.read_visium(path=spatial_path, count_file=h5_path, load_images=True)
    adata_vis.var_names_make_unique()
    adata_vis.obs["sample"] = sample
    slides[sample] = adata_vis

adata_vis = sc.concat(list(slides.values()), label="sample", keys=visium_samples)
```

### Cell 9 — H&E plot (Figure 2c) ✅
```python
sc.pl.spatial(slides["ST8059048"], img_key="hires", color=None,
              title="H&E - ST8059048", save="figure_2c_HE.png")
```

### Cell 10 — Prepare snRNA-seq reference, filter genes ✅
```python
adata_ref = adata_snrna.copy()
adata_ref.X = adata_ref.layers["counts"].copy()

# Manual gene filtering — keep genes expressed in >5% of cells
n_cells = adata_ref.shape[0]
gene_counts = (adata_ref.X > 0).sum(axis=0)
if hasattr(gene_counts, 'A1'):
    gene_counts = gene_counts.A1
pct_expressed = gene_counts / n_cells
selected = pct_expressed > 0.05
adata_ref = adata_ref[:, selected].copy()
```

### Cell 11 — Estimate reference signatures via run_regression ⏳ RUNNING (~5 hours on CPU)
```python
adata_ref.obs['dummy_covar'] = '1'

results_regression = cell2location.run_regression(
    sc_data=adata_ref,
    train_args={
        'covariate_col_names': ['dummy_covar'],
        'sample_name_col': 'sample',
        'tech_name_col': None,
        'n_epochs': 250,
        'minibatch_size': 2500,
        'learning_rate': 0.01,
        'use_average_as_initial_value': True,
        'use_cuda': False,
        'train_proportion': 0.9,
        'l2_weight': True,
        'use_raw': False,
    },
    export_args={'path': './results_regression/', 'save_model': False}
)

inf_aver = results_regression['mod'].samples['post_sample_means']['per_cluster_mu_fg']
inf_aver = pd.DataFrame(
    inf_aver,
    index=adata_ref.var_names,
    columns=adata_ref.uns['mod_cell_types']
)
```

### Cells 12–15 — NOT YET RUN
Pending completion of Cell 11. Plan:
- Cell 12: Find shared genes between Visium and reference signatures
- Cell 13: Run `cell2location.run_cell2location` with `N_cells_per_location=8`,
  `detection_alpha=200`, 20,000 epochs
- Cell 14: Export posterior abundances
- Cell 15: Plot Figures 2d, 2e, 2f

---

## Difficulties and How We Overcame Them

### 1. Job scheduler confusion — sbatch vs qsub
**Problem:** Initially assumed BU's SCC2 used SLURM (`sbatch`) for job submission,
consistent with many other HPC clusters. The `sbatch` command was not available and
failed silently or with errors.
**Solution:** BU's SCC2 uses the SGE/PBS scheduler, where the correct job submission
command is `qsub`. However, given the Singularity constraint (see below), job
submission was ultimately not used — everything runs interactively.

### 2. Singularity only works in desktop interactive session
**Problem:** `singularity shell` failed in regular terminal login sessions.
**Solution:** Must open a desktop interactive session through the HPC web portal first,
then run `singularity shell cell2location.sif` inside that session's terminal, then
launch `jupyter notebook` from inside the container.

### 3. Attempted separate .py scripts — NameError
**Problem:** Tried running a separate `figure2b.py` script, got
`NameError: name 'adata_snrna' is not defined` because the script has no access to
notebook memory.
**Solution:** Abandoned all separate scripts. Everything is in one continuous notebook.

### 4. Old cell2location API in container
**Problem:** The singularity container uses an older version of cell2location that has
no `cell2location.utils`, no `RegressionModel` class, and no `Cell2location` class.
Attempts to import these raised `ModuleNotFoundError` and `ImportError`.
**Solution:** Switched entirely to the function-based old API:
`cell2location.run_regression()` and `cell2location.run_cell2location()`. Identified
the correct API by running `dir(cell2location)` to inspect available functions.

### 5. filter_genes not available
**Problem:** `from cell2location.utils.filtering import filter_genes` failed because
`cell2location.utils` does not exist in this version.
**Solution:** Implemented manual gene filtering using numpy directly — kept genes
expressed in >5% of cells, which approximates the paper's default filtering behavior.

### 6. Double spatial path nesting
**Problem:** `sc.read_visium` raised `OSError: Could not find tissue_positions_list.csv`
because the path resolved to `ST8059048_spatial/spatial/spatial/` — the function
appends `spatial/` internally.
**Solution:** Changed `spatial_path` to point to `ST8059048_spatial/` instead of
`ST8059048_spatial/spatial/`, letting `sc.read_visium` handle the final `spatial/`
append itself.

### 7. squidpy import error
**Problem:** Initial Cell 8 included `import squidpy` which failed because squidpy
is not in the conda environment or container.
**Solution:** Removed squidpy entirely — it was unnecessary. `sc.read_visium` and
`sc.pl.spatial` are pure scanpy functions.

### 8. run_regression AttributeError on adata.raw
**Problem:** `run_regression` raised `AttributeError: 'NoneType' object has no
attribute 'X'` because it defaulted to `use_raw=True` but `adata_ref.raw` was None.
**Solution:** Added `'use_raw': False` to `train_args` to force use of `adata_ref.X`
directly, which contains the raw counts we set explicitly.

### 9. covariate_col_names cannot be [None]
**Problem:** Passing `'covariate_col_names': [None]` caused a crash inside
`run_regression` when it tried to access that column in `adata_ref.obs`.
**Solution:** Added a dummy constant column to `adata_ref.obs`:
`adata_ref.obs['dummy_covar'] = '1'` and passed `'covariate_col_names': ['dummy_covar']`.
This satisfies the library's requirement for a covariate column without affecting the
model since all cells have the same value.

### 10. No GPU available
**Problem:** The desktop interactive session does not have a GPU attached, so
`use_cuda=True` would fail and all training runs on CPU only.
**Impact:** `run_regression` with 250 epochs is estimated to take ~5 hours. The
subsequent `run_cell2location` with 20,000 epochs will take significantly longer.
**Mitigation:** Set `use_cuda=False` explicitly. Consider requesting a GPU-enabled
interactive session from the cluster administrators via `qsub`, or reducing epochs
for a test run first.

---

## Runtime Notes
- Cell 6 (UMAP): moderate runtime, completed successfully
- Cell 11 (run_regression, 250 epochs, CPU): estimated ~5 hours
- Cell 13 (run_cell2location, 20,000 epochs, CPU): estimated many hours — plan ahead

---

## Key Technical Decisions
| Decision | Reason |
|----------|--------|
| Use `layers["counts"]` for cell2location input | cell2location requires raw unnormalized integer counts; `adata.X` was overwritten by normalization in Cell 6 |
| Use `annotation_1` not `annotation_1_print` | Clean label without number prefix; `annotation_1_print` is for figure ordering only |
| BBKNN with `neighbors_within_batch=3` | Matches paper methods exactly |
| `N_cells_per_location=8` | From paper Table 1 for mouse brain Visium |
| `detection_alpha=200` | From paper Table 1 for mouse brain Visium |
| 20,000 epochs for spatial mapping | From paper Table 1 for mouse brain Visium |
| dummy_covar column | Required by old API — constant value has no effect on model |
| `use_raw=False` | `adata_ref.raw` is None; counts are in `adata_ref.X` directly |

---

## Files Produced
| File | Contents | Status |
|------|----------|--------|
| `figure_2b.png` | UMAP with 59 cell subtypes and broad class | ✅ Done |
| `figure_2c_HE.png` | H&E histology image of representative Visium section | ✅ Done |
| `results_regression/` | Output directory from run_regression | ⏳ Running |
| `figure_2d.png` | Spatial map of major regional subtypes | ⏳ Pending |
| `figure_2e.png` | Spatial map of sparse inhibitory neurons | ⏳ Pending |
| `figure_2f.png` | Spatial map of cortical excitatory layer subtypes | ⏳ Pending |
