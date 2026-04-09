# Cell2location Analysis of Human Glioblastoma: Complete Paper Trail

**Project:** DS 596 Special Topics: Learning from Large-scale Biological Data  
**Team:** Team 7 (Vaidehi Gupta, Chloe Ploss, Addison Yam)  
**Last Updated:** 2026-04-09  
**Objective:** Apply Cell2location to map cell types in human glioblastoma Visium data using GBmap as a single-cell reference atlas

---

## 1. Data Acquisition

### 1.1 Reference Dataset: GBmap (GSE211376)

**Source:** Ruiz-Moreno et al. 2022 - Harmonized atlas of 1.1 million GBM cells  
**GEO Accession:** [GSE211376](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211376)  
**Download Location:** `/projectnb/ds596/projects/Team_7/data/GBmap_reference/`

**Download Commands:**
```bash
cd /projectnb/ds596/projects/Team_7/data
mkdir GBmap_reference
cd GBmap_reference

# Download raw counts matrix (genes x cells)
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211376/suppl/GSE211376_raw_counts_Ruiz2022_all_samples_filtered_cells.tsv.gz

# Download cell metadata with annotations
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211376/suppl/GSE211376_metadata_Ruiz2022_all_samples_filtered_cells.csv.gz
```

| File | Size | Description |
|------|------|-------------|
| GBmap_raw_counts.tsv.gz | ~500 MB | Raw expression counts (27,102 genes × 39,355 cells) |
| GBmap_metadata.csv.gz | ~10 MB | Cell metadata with `predicted.high_hierarchy` annotations |

---

### 1.2 Spatial Dataset: Human Glioblastoma Visium

**Source:** 10X Genomics - Parent Visium Human Glioblastoma dataset  
**Download Location:** `/projectnb/ds596/projects/Team_7/data/humanGlioblastoma/`

**Download Commands:**
```bash
cd /projectnb/ds596/projects/Team_7/data
mkdir humanGlioblastoma
cd humanGlioblastoma

# Download filtered count matrix
wget https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5

# Download spatial folder with tissue images and positions
wget https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma_spatial.tar.gz
tar -xzf Parent_Visium_Human_Glioblastoma_spatial.tar.gz
```

| File | Size | Description |
|------|------|-------------|
| Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5 | ~50 MB | Filtered count matrix |
| spatial/ directory | ~10 MB | Tissue images (hires/lowres) and spot coordinates |

---

## 2. Environment Setup

### 2.1 SCC GPU Session Configuration

**Resource Request (SCC OnDemand Desktop):**

| Parameter | Value |
|-----------|-------|
| Partition | l40s (L40S GPU) |
| GPUs | `-l gpus=1` |
| Cores | 4 |
| Memory | 64 GB |
| Time limit | 4 hours |

**GPU Available:**

| GPU | Memory | Status |
|-----|--------|--------|
| L40S | 48 GB | Used for training |

---

### 2.2 Container Setup

```bash
cd /projectnb/ds596/projects/Team_7
singularity shell --nv cell2location.sif
```

**GPU Verification:**
```bash
python -c "import torch; print('CUDA available:', torch.cuda.is_available())"
# Output: CUDA available: True
```

---

### 2.3 Required Python Packages

```
scanpy>=1.9.0
anndata>=0.8.0
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.4.0
cell2location>=0.1.0  (container version)
torch>=1.10.0         (with CUDA support)
```

---

## 3. Data Loading and Preprocessing

### 3.1 Load GBmap Reference

```python
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

ref_dir = "/projectnb/ds596/projects/Team_7/data/GBmap_reference"

# Load metadata
metadata = pd.read_csv(
    os.path.join(ref_dir, "GBmap_metadata.csv.gz"),
    compression='gzip'
)

# Load counts
counts = pd.read_csv(
    os.path.join(ref_dir, "GBmap_raw_counts.tsv.gz"),
    compression='gzip', sep='\t', index_col=0
)

# Create AnnData (cells x genes)
adata_ref = ad.AnnData(
    X=counts.T.values.astype(np.int32),
    obs=metadata,
    var=pd.DataFrame(index=counts.index)
)

# Add clean cell type column
adata_ref.obs['cell_type'] = adata_ref.obs['predicted.high_hierarchy']

# Basic filtering
sc.pp.filter_cells(adata_ref, min_genes=200)
sc.pp.filter_genes(adata_ref, min_cells=3)
```

> **Result:** 39,355 cells × 27,093 genes → filtered to 39,355 cells × 25,731 genes

---

### 3.2 Load Visium Spatial Data

```python
visium_dir = "/projectnb/ds596/projects/Team_7/data/humanGlioblastoma"

adata_spatial = sc.read_visium(
    visium_dir,
    count_file='Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5',
    load_images=True
)
```

> **Result:** 3,468 spots × 36,601 genes

---

### 3.3 Align Reference and Spatial Data

```python
# Make gene names unique (handles duplicate gene names)
adata_ref.var_names_make_unique()
adata_spatial.var_names_make_unique()

# Find common genes
common_genes = np.intersect1d(adata_spatial.var_names, adata_ref.var_names)
print(f"Common genes: {len(common_genes)}")
# Output: 1,223 common genes

# Subset both to common genes
adata_ref = adata_ref[:, common_genes].copy()
adata_spatial = adata_spatial[:, common_genes].copy()

# Save raw counts for Cell2location
adata_spatial.layers["counts"] = adata_spatial.X.copy()
```

---

## 4. Reference Signature Creation

### 4.1 Cell Type Distribution

| Cell Type | Count | Percentage |
|-----------|-------|------------|
| Oligodendrocyte | 12,716 | 32.3% |
| OPC-like | 5,201 | 13.2% |
| TAM-MG | 4,011 | 10.2% |
| TAM-BDM | 3,215 | 8.2% |
| MES-like | 2,790 | 7.1% |
| AC-like | 2,752 | 7.0% |
| NPC-like | 2,687 | 6.8% |
| Neuron | 2,325 | 5.9% |
| OPC | 1,587 | 4.0% |
| Astrocyte | 1,091 | 2.8% |
| Endothelial | 394 | 1.0% |
| CD4/CD8 | 353 | 0.9% |
| Pericyte | 224 | 0.6% |
| DC | 9 | 0.02% |

---

### 4.2 Create Reference Signatures (Mean Expression per Cell Type)

```python
# Ensure cell_type is categorical
adata_ref.obs['cell_type'] = adata_ref.obs['cell_type'].astype('category')
cell_types = adata_ref.obs['cell_type'].cat.categories.tolist()

# Normalize for averaging
sc.pp.normalize_total(adata_ref, target_sum=1e4)
sc.pp.log1p(adata_ref)

# Calculate mean expression per cell type
mean_expression = []
for ct in cell_types:
    cells_in_type = adata_ref[adata_ref.obs['cell_type'] == ct]
    mean_exp = cells_in_type.X.mean(axis=0)
    if hasattr(mean_exp, 'A1'):
        mean_exp = mean_exp.A1
    mean_expression.append(mean_exp)

# Create reference signatures dataframe
ref_signatures = pd.DataFrame(
    np.array(mean_expression).T,
    index=adata_ref.var_names,
    columns=cell_types
)

# Save
ref_signatures.to_csv(os.path.join(ref_dir, "reference_signatures.csv"))
```

> **Result:** 1,223 genes × 14 cell types

---

## 5. Spatial Mapping with Cell2location

### 5.1 Model Configuration

```python
from cell2location.models import LocationModelLinearDependentW

# Prepare data
X_spatial = adata_spatial.X.copy()
if hasattr(X_spatial, 'toarray'):
    X_spatial = X_spatial.toarray()

cell_state_mat = ref_signatures_aligned.values  # Shape: (1223, 14)

# Create model
mod_spatial = LocationModelLinearDependentW(
    cell_state_mat=cell_state_mat,
    X_data=X_spatial,
    n_comb=50,
    data_type='float32',
    n_iter=1000,                    # Training iterations
    learning_rate=0.01,             # Learning rate
    total_grad_norm_constraint=200,
    verbose=True,
    cell_number_prior={
        'cells_per_spot': 15,
        'factors_per_spot': 10,
        'combs_per_spot': 3
    }
)

mod_spatial.fit_advi_iterative()
```

---

### 5.2 Training

**Training Log:**
```
Training spatial mapping model...
Total iterations: 1000
Estimated time: 2-4 hours
100.00% [1000/1000 1:10:46<00:00 Average Loss = 6.0297e+06]
Finished [100%]: Average Loss = 6.0097e+06
```

**Training Summary:**

| Metric | Value |
|--------|-------|
| Iterations completed | 1,000/1,000 |
| Total training time | 1 hour 10 minutes |
| Final loss | 6.01e+06 |
| GPU used | Yes (L40S) |

---

### 5.3 Posterior Sampling (Partially Completed)

```python
mod_spatial.sample_posterior()
```

**Sampling Log:**
```
79.50% [795/1000 56:10<14:29 Average Loss = 6.606e+06]
```

> **Sampling Status:** Interrupted at 79.5% due to GPU idle timeout (2 hours)

---

## 6. Technical Challenges and Solutions

### 6.1 GPU Idle Timeout Issue

**Problem:** SCC terminated job because GPU remained idle for 2 hours during CPU-only posterior sampling.  
**Root Cause:** The `sample_posterior()` method in this Cell2location version runs on CPU only.

**Proposed Solutions:**
1. Split into two jobs: GPU for training, CPU for sampling
2. Reduce `num_samples` to 100–200
3. Skip sampling (use mean estimates only)

---

### 6.2 Cell2location API Version Differences

**Issue:** Multiple API incompatibilities with the container version.

**Resolved Issues:**
- `filter_genes` not available → implemented manual filtering
- `setup_anndata` method missing → used direct model initialization
- `train()` method missing → used `fit_advi_iterative()`
- `export2adata` method missing → manual extraction planned

---

### 6.3 Data Alignment

**Issue:** Variable names not unique in Visium data.  
**Solution:** Used `var_names_make_unique()` before subsetting.

---

## 7. Output Files Generated

| File | Location | Size | Status |
|------|----------|------|--------|
| GBmap_raw_counts.tsv.gz | data/GBmap_reference/ | 500 MB | ✅ |
| GBmap_metadata.csv.gz | data/GBmap_reference/ | 10 MB | ✅ |
| reference_signatures.csv | data/GBmap_reference/ | 15 MB | ✅ |
| spatial_mapping_results.h5ad | results/ | Pending | ⏳ |
| cell2location_results.h5ad | results/ | Pending | ⏳ |

---

## 8. Next Steps

### 8.1 Immediate Actions

1. **Resume from interrupted sampling (reduced samples):**
```python
mod_spatial.sample_posterior(num_samples=100)
```

2. **Extract results:**
```python
# Extract abundance estimates
cell_types = ref_signatures_aligned.columns.tolist()
for i, ct in enumerate(cell_types):
    adata_spatial.obs[f'{ct}_mean'] = mod_spatial.samples['post_sample_means']['spot_factors'][:, i]

# Identify dominant cell type
abundance_cols = [f'{ct}_mean' for ct in cell_types]
dominant_idx = np.argmax(adata_spatial.obs[abundance_cols].values, axis=1)
adata_spatial.obs['dominant_cell_type'] = [cell_types[i] for i in dominant_idx]
```

---

### 8.2 Figure Generation (Planned)

| Figure | Description | Status |
|--------|-------------|--------|
| 2b | UMAP of reference cell types | Pending |
| 2c | Cell type composition bar plot | Pending |
| 2d–f | Spatial mapping of GBM cell states | Pending |
| Extended | Dominant cell type map | Pending |
| Extended | Regional heterogeneity analysis | Pending |

---

### 8.3 Expected Output Directory Structure

```
/projectnb/ds596/projects/Team_7/
├── data/
│   ├── GBmap_reference/
│   │   ├── GBmap_raw_counts.tsv.gz
│   │   ├── GBmap_metadata.csv.gz
│   │   └── reference_signatures.csv
│   └── humanGlioblastoma/
│       ├── Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5
│       └── spatial/
├── results/
│   ├── spatial_mapping_results.h5ad
│   ├── cell2location_results.h5ad
│   ├── cell_abundance_summary.csv
│   ├── spot_abundance_matrix.csv
│   └── spot_annotations.csv
├── figures/
│   ├── fig2b_reference_umap.png
│   ├── fig2c_composition.png
│   ├── fig2def_spatial_mapping.png
│   ├── dominant_cell_type_map.png
│   └── figure2_complete.png
└── notebooks/
    ├── 01_data_download.ipynb
    ├── 02_reference_preparation.ipynb
    └── 03_cell2location_analysis.ipynb
```

---

## 9. References

1. Kleshchevnikov, V., Shmatko, A., Dann, E. et al. Cell2location maps fine-grained cell types in spatial transcriptomics. *Nat Biotechnol* 40, 661–671 (2022).
2. Ruiz-Moreno, et al. Harmonized single-cell atlas of glioblastoma. GSE211376 (2022).
3. Neftel, C., Laffy, J., Filbin, M.G. et al. An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. *Cell* 178(4):835–849.e21 (2019).
4. 10X Genomics. Parent Visium Human Glioblastoma dataset. (2020).

---

## 10. Session Information

```
Date:           2026-04-09
Container:      cell2location.sif
GPU:            NVIDIA L40S (48GB)
Python:         3.7.8
PyTorch:        with CUDA support
Scanpy:         1.9.0
Cell2location:  version in container
```

---

*End of Paper Trail*
