# Cell2Location: Reproduction and Extensions to Glioblastoma and Xenium

**DS 596 Special Topics: Learning from Large-scale Biological Data**  
**Team 7** — Vaidehi Gupta, Chloe Ploss, Addison Yam  
**Professor:** Pawel Przytycki  
**Based on:** Kleshchevnikov et al. (2022), *Nature Biotechnology*  
**Full Report:** DS596_Final_Report.pdf  
**GitHub:** https://github.com/cploss-dot/Cell2Location-Project-Research-Reproduction-Extension

---

## Table of Contents

- [Project Overview](#project-overview)
- [Repository Structure](#repository-structure)
- [Datasets](#datasets)
  - [Part 1: Reproduction — Mouse Brain Visium](#part-1-reproduction--mouse-brain-visium)
  - [Part 2: Extension — Human Glioblastoma Visium](#part-2-extension--human-glioblastoma-visium)
  - [Part 3: Extension — Xenium Mouse Brain](#part-3-extension--xenium-mouse-brain)
- [Computational Environment](#computational-environment)
- [Environment Setup Guide](#environment-setup-guide)
- [Running the Notebooks](#running-the-notebooks)
- [Notebook Descriptions](#notebook-descriptions)
  - [reference.ipynb](#referenceipynb)
  - [cell2location_reproduction.ipynb](#cell2location_reproductionipynb)
  - [cell2location_visium_human_gbm.ipynb](#cell2location_visium_human_gbmipynb)
  - [cell2location_xenium_mouse_brain.ipynb](#cell2location_xenium_mouse_brainipynb)
- [Results Summary](#results-summary)
- [Known Issues and Limitations](#known-issues-and-limitations)
- [Team Contributions](#team-contributions)
- [References](#references)

---

## Project Overview

This project reproduces key findings from the Cell2location paper (Kleshchevnikov et al., 2022) and extends the method to two new contexts: mapping cell types in human glioblastoma and applying the model to high-resolution single-cell spatial data from 10x Xenium. All work was performed using cell2location v0.3 (PyMC3/Theano-based) inside a Singularity container on Boston University's Shared Computing Cluster (SCC).

**Original Method:** Cell2location is a Bayesian hierarchical model that maps cell types onto spatial transcriptomic data. Step 1 fits a negative binomial regression to a scRNA-seq reference to learn cell type-specific gene expression signatures. Step 2 uses those signatures to deconvolve each spatial spot into estimated cell type abundances, borrowing statistical strength across spots via a hierarchical factorization prior.

**Reproductions:**  
Figures 2b, 2c, 2d, and 2e from the original paper were successfully reproduced using the authors' publicly available mouse brain datasets. Two training durations (3,000 and 30,000 iterations) were compared to demonstrate the effect of model convergence on spatial resolution.

**Extensions:**
- **Human Glioblastoma Visium:** Cell2location applied to human glioblastoma Visium data using the GBmap harmonized atlas (39,355 cells, 14 cell types) as reference. Malignant states (MES-like, NPC-like, OPC-like, AC-like) and tumor microenvironment populations were spatially mapped. A parameter sensitivity analysis was performed across five configurations.
- **Xenium Mouse Brain:** Cell2location applied to 10x Xenium data at near single-cell resolution (~10 µm per location). Five parameter configurations were tested to characterize how the resolution change affects model outputs and what adjustments are necessary.

---

## Repository Structure

```
Cell2Location-Project-Research-Reproduction-Extension/
│
├── data/
│   ├── cell2location.sif               ← Singularity container image (~4.8GB)
│   ├── mousescRNAseq/                  ← Mouse brain snRNA-seq (E-MTAB-11115)
│   ├── mouseVisium/                    ← Mouse brain Visium (E-MTAB-11114)
│   ├── humanGlioblastoma/              ← Human GBM Visium spatial data
│   ├── GBmap_reference/                ← Processed GBmap scRNA-seq reference
│   └── xeniumMouseBrain/               ← Xenium mouse brain coronal subset
│
├── notebooks/
│   ├── cell2location_reproduction.ipynb
│   ├── cell2location_xenium_mouse_brain.ipynb
│   ├── cell2location_visium_human_gbm.ipynb
│   └── reference.ipynb
│
├── reproduction_figures/               ← All reproduction output figures
├── human_gbm_extension_results/        ← All GBM extension output figures
├── xenium_extension_results/           ← All Xenium extension output figures
├── project_environment.yml             ← Conda environment record (not used at runtime)
├── DS596_Final_Report.pdf
└── README.md
```

---

## Datasets

### Part 1: Reproduction — Mouse Brain Visium

#### snRNA-seq Reference (E-MTAB-11115)

Six adult mouse brain snRNA-seq samples. 40,531 annotated cells across 59 fine-grained cell types. Downloaded from EBI BioStudies:

```bash
wget -r --no-parent \
  ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/115/E-MTAB-11115/Files/
mv ftp.ebi.ac.uk/.../Files/ mousescRNAseq/
```

**Files used:**

| File | Description |
|---|---|
| `5705STDY8058280_filtered_feature_bc_matrix.h5` | Sample 1 |
| `5705STDY8058281_filtered_feature_bc_matrix.h5` | Sample 2 |
| `5705STDY8058282_filtered_feature_bc_matrix.h5` | Sample 3 |
| `5705STDY8058283_filtered_feature_bc_matrix.h5` | Sample 4 |
| `5705STDY8058284_filtered_feature_bc_matrix.h5` | Sample 5 |
| `5705STDY8058285_filtered_feature_bc_matrix.h5` | Sample 6 |
| `cell_annotation.csv` | 40,531 cell IDs mapped to `annotation_1` labels |

**cell_annotation.csv format:**
```
Cell ID,sample,annotation_1,annotation_1_print
5705STDY8058283_AAACCCAAGCCTATTG-1,5705STDY8058283,Ext_L23,22_Ext_L23
```
Cell ID format: `{sampleID}_{barcode}-1`

**Files not used:** `*_raw_feature_bc_matrix.h5`, `*_molecule_info.h5`, `*_web_summary.html`, `*_metrics_summary.csv`

#### Visium Spatial Data (E-MTAB-11114)

Five adult mouse brain coronal sections across two animals. Spatial mapping was performed on section ST8059048 (2,987 spots, 10,085 shared genes with reference).

```bash
wget -r --no-parent \
  ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/114/E-MTAB-11114/Files/
mv ftp.ebi.ac.uk/.../Files/ mouseVisium/
```

**Post-download — extract all spatial archives:**

```bash
cd mouseVisium/
for f in *.tar.gz; do
    sample="${f%%_spatial.tar.gz}"
    mkdir -p "${sample}_spatial"
    tar -xzf "$f" -C "${sample}_spatial"
done
```

> **Important:** Each spatial path has an extra nesting level — it is `ST8059048_spatial/spatial/` not `ST8059048_spatial/`. The `read_visium` call in the notebook points to the inner `spatial/` directory.

Each extracted `spatial/` folder contains:
```
aligned_fiducials.jpg
detected_tissue_image.jpg
scalefactors_json.json
tissue_hires_image.png
tissue_lowres_image.png
tissue_positions_list.csv
```

**Files used for spatial mapping:**
- `ST8059048_filtered_feature_bc_matrix.h5`
- `ST8059048_spatial/spatial/`

---

### Part 2: Extension — Human Glioblastoma Visium

#### scRNA-seq Reference (GBmap — GSE211376)

GBmap harmonized single-cell atlas: 39,355 cells from 110 glioblastoma patients, 14 annotated cell populations. Downloaded and processed by `reference.ipynb` into `GBmap_reference/`.

```bash
mkdir -p data/GBmap_reference && cd data/GBmap_reference

wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211376/suppl/GSE211376_raw_counts_Ruiz2022_all_samples_filtered_cells.tsv.gz" \
    -O GBmap_raw_counts.tsv.gz

wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211376/suppl/GSE211376_metadata_Ruiz2022_all_samples_filtered_cells.csv.gz" \
    -O GBmap_metadata.csv.gz
```

Cell type labels come from the `predicted.high_hierarchy` column in the metadata. After QC filtering (min 200 genes/cell, min 3 cells/gene) and gene intersection with Visium data: **25,731 shared genes** retained.

**Cell type distribution:**

| Cell Type | Count | % |
|---|---|---|
| Oligodendrocyte | 12,716 | 32.3% |
| OPC-like | 5,201 | 13.2% |
| TAM-MG | 4,011 | 10.2% |
| TAM-BDM | 3,215 | 8.2% |
| MES-like | 2,790 | 7.1% |
| AC-like | 2,752 | 7.0% |
| NPC-like | 2,687 | 6.8% |
| + 7 more types | ~5,983 | ~15.2% |

**Processed output saved to:**
- `GBmap_reference/GBmap_reference_compressed.h5ad` — loaded directly by the GBM notebook on subsequent runs

#### Visium Spatial Data (Human Glioblastoma)

10x Genomics human glioblastoma Visium: 3,468 spots (all in tissue), 36,601 genes.

```bash
mkdir data/humanGlioblastoma && cd data/humanGlioblastoma

wget https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5

wget https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma_spatial.tar.gz

tar -xzf Parent_Visium_Human_Glioblastoma_spatial.tar.gz
```

**Quick sanity check:**
```python
import scanpy as sc, numpy as np
adata = sc.read_visium(
    path="data/humanGlioblastoma/",
    count_file="Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5"
)
adata.var_names_make_unique()
print(adata.shape)                                              # (3468, 36601)
print(np.all(adata.X.data == adata.X.data.astype(int)))        # True
print(adata.obs["in_tissue"].sum())                             # 3468
```

**Parameter sweep configurations:**

| Run | N_iter_reg | N_iter_spat | N_comb | LR | N_samples | Purpose |
|---|---|---|---|---|---|---|
| 1 Test | 100 | 500 | 20 | 0.005 | 200 | Baseline |
| 2 More Samples | 100 | 500 | 20 | 0.005 | 500 | Sampling effects |
| 3 Higher LR | 100 | 500 | 20 | 0.025 | 200 | Convergence speed |
| 4 Intermediate | 500 | 5,000 | 30 | 0.005 | 500 | Balanced |
| 5 Final | 2,000 | 10,000 | 30 | 0.005 | 500 | Paper quality |

Cell number prior: **15 cells/spot** (derived from histology image inspection, adjusted for GBM's high cellularity).

---

### Part 3: Extension — Xenium Mouse Brain

#### scRNA-seq Reference (Cell2location Tutorial — Mouse Brain)

Official cell2location mouse brain reference atlas: 40,572 annotated cells, full transcriptome (31,053 genes). **Downloaded automatically by the Xenium notebook** from the Bayraktar Lab's public hosting:

```python
# Auto-downloaded to RESULTS_DIR if not already present:
# https://cell2location.cog.sanger.ac.uk/tutorial/mouse_brain_snrna/all_cells_20200625.h5ad
# https://cell2location.cog.sanger.ac.uk/tutorial/mouse_brain_snrna/snRNA_annotation_astro_subtypes_refined59_20200823.csv
```

After converting reference gene names from Ensembl IDs to gene symbols and intersecting with the Xenium panel: **247 shared genes** across 59 cell types. After QC filtering (min 10 counts/cell, min 5 cells/gene): 183 cells removed, 0 genes removed.

#### Xenium Spatial Data (Mouse Brain Coronal Subset)

Near single-cell resolution (~10 µm) spatial transcriptomics of mouse brain cortex (CTX) and hippocampus (HP). 36,602 cells, 248 genes, mean 246.7 transcripts/cell.

```bash
mkdir data/xeniumMouseBrain && cd data/xeniumMouseBrain

wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip

unzip Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
```

Dataset page: https://www.10xgenomics.com/datasets/xenium-ff-mouse-brain-coronal-subset-ctx-hp-1-0-2

**Key files used:**

| File | Description |
|---|---|
| `cell_feature_matrix.h5` | Cell × gene count matrix (36,602 × 248) |
| `cells.csv.gz` | x/y centroids, transcript counts, cell area, nucleus area |
| `morphology.ome.tif` | High-resolution fluorescence tissue image |
| `gene_panel.json` | 275 targeted genes (248 after QC) |

**Parameter sweep configurations:**

| | | Run1 | Run2 | Run3 | Run4 | Run5 |
|---|---|---|---|---|---|---|
| **Train Model** | N_iter_reg | 100 | 500 | 100 | 100 | 4,000 |
| | Learning_rate | 0.005 | 0.005 | 0.025 | 0.005 | 0.005 |
| | Train/Test Split | 90/10% | 90/10% | 90/10% | 70/30% | 90/10% |
| **Spatial Model** | N_iter_spatial | 500 | 2,500 | 500 | 500 | 20,000 |
| | N_comb | 20 | 20 | 20 | 20 | 50 |
| | Learning_rate | 0.005 | 0.005 | 0.025 | 0.005 | 0.005 |
| **Sample Posterior** | N_samples | 200 | 200 | 200 | 400 | 200 |
| | Batch_size | 200 | 200 | 200 | 400 | 200 |

---

## Computational Environment

All analyses were performed on Boston University's Shared Computing Cluster (SCC2).

**Critical constraints:**

| Constraint | Detail |
|---|---|
| Cell2location version | v0.3 — PyMC3/Theano backend. Current v0.1.3 uses PyTorch/scvi-tools with a different API |
| Why v0.3 | Newer version dependencies exceed 10GB SCC project storage quota |
| GPU required | Theano backend requires GPU. **Only compatible with V100 (Volta) — incompatible with A100 and L40s** |
| GPU request order | GPU must be requested before launching Singularity — cannot be added after the container starts |
| `sbatch` unavailable | No job scheduler — everything runs interactively |
| Singularity restriction | Container only works in a desktop interactive session, not a regular login terminal |
| Theano GPU flag | Must set `THEANO_FLAGS = "device=cuda,floatX=float32"` as the absolute first cell in every notebook, before any imports |

**Working directory on SCC:**
```
/projectnb/ds596/projects/Team_7/data/
```

---

## Environment Setup Guide

### Pull the Container (Run Once)

```bash
cd /projectnb/ds596/projects/Team_7/data/
singularity pull cell2location.sif docker://quay.io/vitkl/cell2location:latest
```

Creates a single ~4.8GB `.sif` file.

### Launch Workflow (Every Session)

**Step 1 — Request a Desktop interactive session via SCC OnDemand:**
- App: Desktop
- GPU: Request **Volta/V100 architecture specifically**
- RAM: At least 16GB
- Time: Reproduction ~2h, GBM full run ~4h, Xenium Run5 ~3h

**Step 2 — Open a terminal in the desktop and verify GPU:**
```bash
nvidia-smi
# Must show a V100
```

**Step 3 — Launch the container with GPU passthrough:**
```bash
singularity exec --nv /projectnb/ds596/projects/Team_7/data/cell2location.sif bash
```
The `--nv` flag is mandatory. Without it the GPU is inaccessible.

**Step 4 — Start Jupyter from inside the container:**
```bash
jupyter lab --ip=0.0.0.0 --port=8888 --no-browser
# or: jupyter notebook
```

**Step 5 — Connect:**
Copy the URL with token from the terminal (e.g. `http://127.0.0.1:8888/lab?token=...`) and open it in the desktop session browser.

**Step 6 — MANDATORY first cell in every notebook:**
```python
import os
os.environ["THEANO_FLAGS"] = "device=cuda,floatX=float32"
```
This must execute before any `import cell2location`, `import pymc3`, or `import theano`. Skipping it causes the model to run on CPU, making training prohibitively slow.

### Conda Environment (Not Used — For Reference Only)

```bash
conda env create -f project_environment.yml
```
Failed — scvi-tools dependencies for v0.1.3 exceed the 10GB quota. Included for reference only.

---

## Running the Notebooks

**All notebooks must run inside the Singularity container via the desktop session workflow above.** Do not run them as external `.py` scripts — variables must persist between cells.

**Run order:**

| Order | Notebook | Prerequisite | Approx. Runtime (V100 GPU) |
|---|---|---|---|
| 1 | `reference.ipynb` | `GBmap_reference/` folder created, wget commands run | ~15 min |
| 2 | `cell2location_reproduction.ipynb` | `mousescRNAseq/` and `mouseVisium/` downloaded and extracted | ~25 min (30k iter) |
| 3 | `cell2location_visium_human_gbm.ipynb` | `reference.ipynb` completed | ~30 min per configuration |
| 4 | `cell2location_xenium_mouse_brain.ipynb` | `xeniumMouseBrain/` downloaded | ~20 min (Run1) to ~3h (Run5) |

**Session persistence warning:** If the desktop session times out, all in-memory variables are lost and cells must be rerun from the beginning. The GBM notebook saves a compressed reference `.h5ad` after the first load (`GBmap_reference_compressed.h5ad`) to avoid re-loading the large count matrix on subsequent runs.

---

## Notebook Descriptions

### cell2location_reproduction.ipynb

Reproduces Figures 2b through 2e from the original paper.

**Key steps:**

**1. Load snRNA-seq data**
```python
# Loads all 6 samples and concatenates into one AnnData
# Cell IDs formatted as {sampleID}_{barcode} to match cell_annotation.csv
adata_snrna = sc.concat(adatas)
# Merges annotation_1 labels; drops unannotated cells
# Preserves raw counts in adata_snrna.layers["counts"]
```

**2. Figure 2b — UMAP**
```python
# normalize → log1p → HVG (3000 genes, batch_key="sample")
# → scale → PCA (50 components)
# → BBKNN batch correction (neighbors_within_batch=3)
# → UMAP
# Plots: fine subtypes (59) + broad class (7)
```

**3. Figure 2c — H&E**
```python
sc.read_visium(path=spatial_path, count_file=h5_path, load_images=True)
sc.pl.spatial(slides["ST8059048"], img_key="hires", color=None)
```

**4. Reference signature estimation**
```python
# Uses cell2location.run_regression with old PyMC3 API
# Key fix: use dummy_covar column to satisfy covariate_col_names requirement
adata_ref.obs['dummy_covar'] = '1'
results_regression = cell2location.run_regression(
    sc_data=adata_ref,
    train_args={
        'covariate_col_names': ['dummy_covar'],
        'sample_name_col': 'sample',
        'n_epochs': 100,           # test; 250 for full run
        'learning_rate': 0.01,
        'minibatch_size': 2500,
    }
)
inf_aver  # shape: (10,085 genes × 65 cell types)
```

**5. Spatial mapping — Figures 2d and 2e**
```python
# Genes intersected: 10,085 shared between Visium and reference
results_c2l = cell2location.run_cell2location(
    sc_ref=inf_aver,
    sp_data=adata_vis,  # section ST8059048, 2,987 spots
    train_args={
        'n_epochs': 30000,      # 3000 for test run
        'learning_rate': 0.001,
    },
    model_kwargs={
        'cell_number_prior': {
            'cells_per_spot': 8,
            'factors_per_spot': 7,
            'combs_per_spot': 2.5
        },
        'cell_number_var': {'alpha_mean': 200}
    }
)
# Figure 2d: Oligo_2, Inh_Meis2_3, Inh_4, Ext_Thal_1, Ext_L23, Ext_L56
# Figure 2e: Inh_Sst, Inh_Lamp5, Inh_Vip
sc.pl.spatial(adata_plot, color=cell_types, img_key="hires", vmax="p99")
```

**Known API issue:** `covariate_col_names=[None]` and `covariate_col_names=[]` both error in cell2location v0.3. Fix is to pass a dummy constant covariate column.

---

### reference.ipynb

Downloads and processes the GBmap atlas from NCBI GEO for use as the GBM extension reference.

**Key steps:**
```python
# Downloads two files from GSE211376:
# GBmap_raw_counts.tsv.gz   — 27,102 genes × 39,355 cells
# GBmap_metadata.csv.gz     — cell metadata including predicted.high_hierarchy

# Creates AnnData (cells × genes), transposing from gene × cell format
adata_ref = ad.AnnData(X=counts.T.values.astype(np.int32), obs=metadata)
adata_ref.obs['cell_type'] = adata_ref.obs['predicted.high_hierarchy']

# QC filtering
sc.pp.filter_cells(adata_ref, min_genes=200)   # removes low-quality cells
sc.pp.filter_genes(adata_ref, min_cells=3)     # removes sparse genes

# Saves compressed h5ad for fast reloading:
# GBmap_reference/GBmap_reference_compressed.h5ad
```

Run this notebook once before running `cell2location_visium_human_gbm.ipynb`.

---

### cell2location_visium_human_gbm.ipynb

Full GBM extension workflow.

**Key steps:**

**0. GPU setup — must be first cell:**
```python
import os
os.environ["THEANO_FLAGS"] = "device=cuda,floatX=float32"
import theano
print(f"Theano device: {theano.config.device}")
```

**1–4. Load reference and spatial data, align genes:**
```python
# Loads GBmap_reference_compressed.h5ad if exists, else rebuilds from raw files
# Loads Visium with sc.read_visium
# Intersects genes: 25,731 shared genes retained
# Saves raw counts in adata_spatial.layers["counts"]
```

**5. Reference signatures (two methods):**
```python
# Method A — mean expression per cell type (fast)
inf_aver_mean = adata_ref.to_df().groupby(adata_ref.obs['cell_type']).mean().T

# Method B — regression-based (used for spatial mapping)
from cell2location.models import RegressionGeneBackgroundCoverageTorch
mod = RegressionGeneBackgroundCoverageTorch(
    sample_id="patient",
    cell2covar=cell2covar_df,    # sample + cell_type columns
    X_data=X_ref,
    n_iter=N_ITER_REG,           # 100 (test) to 2000 (final)
    learning_rate=0.005,
    use_cuda=True
)
mod.fit_advi_iterative(n=1, n_iter=N_ITER_REG, train_proportion=0.9)
# Extracts inf_aver: (25,731 genes × 14 cell types)
```

**6. Spatial model:**
```python
from cell2location.models import LocationModelLinearDependentW
mod_spatial = LocationModelLinearDependentW(
    cell_state_df=inf_aver,
    X_data=X_spatial,            # 3,468 spots × 25,731 genes
    n_comb=N_COMB,               # 20 (test) to 30 (final)
    n_iter=N_ITER_SPAT,          # 500 (test) to 10,000 (final)
    learning_rate=LR,
    use_cuda=True,
    cell_number_prior={
        'cells_per_spot': 15,    # adjusted for GBM high cellularity
        'factors_per_spot': 7,
        'combs_per_spot': 2.5
    }
)
# Output: spot_factors (3,468 × 14) added to AnnData obs
```

**7. Visualizations per configuration:**
- Dominant cell type spatial map (highest abundance cell type per spot)
- Six key cell types (MES-like, Astrocyte, CD4/CD8, OPC, AC-like, Pericyte)
- All 14 cell types panel

---

### cell2location_xenium_mouse_brain.ipynb

Full Xenium extension workflow.

**Key steps:**

**0. GPU setup — must be first cell:**
```python
import os
os.environ["THEANO_FLAGS"] = "device=cuda,floatX=float32"
```

**1. Load Xenium data:**
```python
adata_xenium = sc.read_10x_h5("xeniumMouseBrain/cell_feature_matrix.h5")
adata_xenium.var_names_make_unique()

cells = pd.read_csv("xeniumMouseBrain/cells.csv.gz").set_index("cell_id")
adata_xenium.obs = adata_xenium.obs.join(cells[["x_centroid", "y_centroid"]])
adata_xenium.obsm["spatial"] = adata_xenium.obs[["x_centroid","y_centroid"]].values
# Shape: (36,602 cells × 248 genes)
```

**2. Download reference (auto if not present):**
```python
# Downloads ~700MB h5ad and annotation CSV from Bayraktar Lab
# Merges annotation_1 cell type labels onto reference
# Swaps Ensembl IDs to gene symbols for gene name matching
adata_ref.var_names = adata_ref.var["SYMBOL"]
```

**3. Find shared genes and QC:**
```python
shared_genes = adata_xenium.var_names.intersection(adata_ref.var_names)
# 247 shared genes
sc.pp.filter_cells(adata_xenium, min_counts=10)   # removes 183 cells
sc.pp.filter_genes(adata_xenium, min_cells=5)     # removes 0 genes
```

**4. Regression model:**
```python
from cell2location.models import RegressionGeneBackgroundCoverageTorch
mod = RegressionGeneBackgroundCoverageTorch(
    sample_id="sample",
    cell2covar=cell2covar,       # sample + annotation_1 columns
    X_data=X_data_ref,
    n_iter=N_ITER_REG,           # 100 (test) to 4000 (Run5)
    learning_rate=0.005,
    use_cuda=False               # set True with GPU
)
mod.fit_advi_iterative(n=1, n_iter=N_ITER_REG, train_proportion=TRAIN_PROPORTION)
# inf_aver: (247 genes × 59 cell types)
```

**5. Spatial model:**
```python
from cell2location.models import LocationModelLinearDependentW
mod_spatial = LocationModelLinearDependentW(
    cell_state_df=inf_aver,
    X_data=X_xenium,             # (36,419 cells × 247 genes) after QC
    n_comb=N_COMB,               # 20 (runs 1-4) to 50 (Run5)
    n_iter=N_ITER_SPAT,          # 500 (Run1) to 20,000 (Run5)
    learning_rate=LR,
    use_cuda=False
)
# Output: spot_factors (36,419 × 59) added to AnnData obs
```

**6. Visualizations:**
```python
# 8 selected cell types (major brain categories):
# Astro_CTX, Oligo_2, Ext_L23, Ext_Hpc_CA1, Inh_Sst, Micro, OPC_1, Endo
# All 59 cell types panel
# Plotted using x_centroid/y_centroid coordinates with viridis colormap
```

**Speed settings (change at top of notebook):**
```python
N_ITER_REG = 100      # test: 100, final: 4000
N_ITER_SPAT = 500     # test: 500, final: 20000
N_COMB = 20           # test: 20,  final: 50
TRAIN_PROPORTION = 0.9
```

---

## Results Summary

### Reproduction

**Figure 2b:** UMAP from 40,532 annotated cells across 65 cell types. Excitatory neurons form the largest continuous manifold. Inhibitory subtypes form distinct peripheral clusters. Topology closely matches the original paper. Minor shape differences expected from stochastic UMAP initialization.

**Figure 2c:** H&E image of section ST8059048 confirmed loaded with correct spatial alignment. Allen Mouse Brain Atlas region overlay omitted (requires tissue registration outside scope of cell2location).

**Figures 2d and 2e:** Run under 3,000 and 30,000 iterations. The 30,000 iteration run shows clear posterior shrinkage — Excitatory L2/3 signal tightens from a diffuse cortex-wide signal to a sharply localized cortical ribbon. Sparse inhibitory neurons (Inh_Sst, Inh_Lamp5, Inh_Vip) show marked spatial concentration at 30,000 iterations versus diffuse signal at 3,000. Supplementary figures S1 and S2 show the 3,000 iteration comparison runs.

### GBM Extension

All five configurations consistently identified MES-like as the dominant malignant state. The test run (100/500 iterations) matched the top-three cell type rankings of the final run (2,000/10,000 iterations), indicating broad spatial patterns are robust to parameter choice. MES-like cells distributed throughout specific tumor regions consistent with their therapy-resistant biology. Endothelial cells correctly localized to vascular structures, providing internal validation. Dominant cell type maps became progressively cleaner over runs, with the final run resolving primarily Astrocyte and MES-like with minor CD4/CD8 and Pericyte contributions.

### Xenium Extension

Run5 (40x more iterations, 2.5x n_comb) produced the smoothest color gradients and closest-to-publication-quality results. Run4 (70/30 train/test split) produced rectangular artifact patches across all 59 cell types — a coordinate/data handling artifact rather than biological signal. Cell2location remained functional at Xenium resolution, with cell types appearing in similar regions and relative abundances as the paper's Visium results for the same brain regions.

---

## Known Issues and Limitations

**Combined cell type figures:** The original paper shows multiple cell types overlaid in one image. Cell2location v0.3 does not include this functionality and the code for it was not found in the authors' GitHub. All figures show each cell type on a separate subplot.

**GPU compatibility:** Theano backend is incompatible with A100 and L40s GPUs. If your cluster no longer supports V100, you must use cell2location v0.1.3 with its entirely different API.

**GPU request order:** The GPU must be requested before starting Singularity with `--nv`. Starting without `--nv` and trying to add the GPU later requires restarting the container.

**Theano flag order:** `os.environ["THEANO_FLAGS"]` must be set before any theano-related imports. If any notebook cell imports these libraries before the flag is set, GPU will not be used and you must restart the kernel.

**covariate_col_names API bug:** In cell2location v0.3, passing `covariate_col_names=[None]` or `covariate_col_names=[]` causes errors. Workaround: add a dummy constant covariate column (`adata_ref.obs['dummy_covar'] = '1'`) and pass `covariate_col_names=['dummy_covar']`.

**Storage constraints:** Newer cell2location v0.1.3 dependencies exceed 10GB. Full pipeline storage (container + all datasets + outputs) requires approximately 20-25GB.

**Session persistence:** Interactive sessions time out. If this happens, all in-memory variables are lost. The GBmap notebook saves a compressed `.h5ad` checkpoint after first load. For reproduction, rerunning from scratch takes approximately 25 minutes.

**Xenium Run4 artifacts:** The rectangular patches of anomalously high abundance in Run4 (70/30 split, 400 posterior samples) persist across all 59 cell types, identifying them as a data handling artifact rather than biology.

**Reference quality dependency:** Cell2location requires high-quality annotated scRNA-seq reference. Ready-to-use references exist only for the three tissues in the original paper. For GBM, a suitable reference had to be found, downloaded, and manually processed — this process significantly impacts result quality and tool adoption.

---

## Team Contributions

**Addison Yam**
- Found and downloaded all datasets
- Produced the human glioblastoma extension including building the GBmap reference from raw GEO files
- Presented methods and GBM extension results
- Wrote the cell2location method overview and GBM extension sections of the report

**Chloe Ploss**
- Solved the GPU configuration workflow within the Singularity container
- Produced the Visium mouse brain reproduction (Figures 2b and 2c) and the Xenium mouse brain extension (all five runs)
- Presented critiques and Xenium extension results
- Wrote the Visium extension methods, results, and discussion sections of the report

**Vaidehi Gupta**
- Selected the primary paper for the project
- Produced the Visium mouse brain reproduction (Figures 2d and 2e)
- Presented the project introduction and reproduction results
- Wrote the introduction, reproduction methods, and results sections of the report

All team members contributed to debugging the SCC, Singularity, and GPU workflow.

---

## References

Kleshchevnikov, V., Shmatko, A., Dann, E. et al. Cell2location maps fine-grained cell types in spatial transcriptomics. *Nat Biotechnol* 40, 661–671 (2022). https://doi.org/10.1038/s41587-021-01139-4

Ruiz-Moreno, C., et al. Harmonized single-cell atlas of glioblastoma. GEO: GSE211376 (2022).

Neftel, C., Laffy, J., Filbin, M.G. et al. An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. *Cell* 178(4):835-849.e21 (2019).

10x Genomics. Parent Visium Human Glioblastoma dataset (2020). https://www.10xgenomics.com/datasets

10x Genomics. Xenium V1 FF Mouse Brain Coronal Subset (CTX + HP) (2023). https://www.10xgenomics.com/datasets/xenium-ff-mouse-brain-coronal-subset-ctx-hp-1-0-2
