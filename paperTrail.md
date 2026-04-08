### Dataset 1: Mouse Brain snRNA-seq (E-MTAB-11115)
```bash
wget -r --no-parent ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/115/E-MTAB-11115/Files/
```
**Output:** Not captured — command was run before logging began.

**Result confirmed by `ls -lh mousescRNAseq/`:**
```
total 2.8G
-rw-r--r-- 1 addisony ds596  27M Nov 17 2023 5705STDY8058280_filtered_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596 683  Nov 17 2023 5705STDY8058280_metrics_summary.csv
-rw-r--r-- 1 addisony ds596 284M Nov 17 2023 5705STDY8058280_molecule_info.h5
-rw-r--r-- 1 addisony ds596 165M Nov 17 2023 5705STDY8058280_raw_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596 4.4M Nov 17 2023 5705STDY8058280_web_summary.html
-rw-r--r-- 1 addisony ds596  28M Nov 17 2023 5705STDY8058281_filtered_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596 683  Nov 17 2023 5705STDY8058281_metrics_summary.csv
-rw-r--r-- 1 addisony ds596 288M Nov 17 2023 5705STDY8058281_molecule_info.h5
-rw-r--r-- 1 addisony ds596 165M Nov 17 2023 5705STDY8058281_raw_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596 4.4M Nov 17 2023 5705STDY8058281_web_summary.html
-rw-r--r-- 1 addisony ds596  22M Nov 17 2023 5705STDY8058282_filtered_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596 683  Nov 17 2023 5705STDY8058282_metrics_summary.csv
-rw-r--r-- 1 addisony ds596 241M Nov 17 2023 5705STDY8058282_molecule_info.h5
-rw-r--r-- 1 addisony ds596 156M Nov 17 2023 5705STDY8058282_raw_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596 4.1M Nov 17 2023 5705STDY8058282_web_summary.html
-rw-r--r-- 1 addisony ds596  22M Nov 17 2023 5705STDY8058283_filtered_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596 683  Nov 17 2023 5705STDY8058283_metrics_summary.csv
-rw-r--r-- 1 addisony ds596 268M Nov 17 2023 5705STDY8058283_molecule_info.h5
-rw-r--r-- 1 addisony ds596 162M Nov 17 2023 5705STDY8058283_raw_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596 4.1M Nov 17 2023 5705STDY8058283_web_summary.html
-rw-r--r-- 1 addisony ds596  15M Nov 17 2023 5705STDY8058284_filtered_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596 683  Nov 17 2023 5705STDY8058284_metrics_summary.csv
-rw-r--r-- 1 addisony ds596 203M Nov 17 2023 5705STDY8058284_molecule_info.h5
-rw-r--r-- 1 addisony ds596 149M Nov 17 2023 5705STDY8058284_raw_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596 3.8M Nov 17 2023 5705STDY8058284_web_summary.html
-rw-r--r-- 1 addisony ds596  34M Nov 17 2023 5705STDY8058285_filtered_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596 684  Nov 17 2023 5705STDY8058285_metrics_summary.csv
-rw-r--r-- 1 addisony ds596 357M Nov 17 2023 5705STDY8058285_molecule_info.h5
-rw-r--r-- 1 addisony ds596 177M Nov 17 2023 5705STDY8058285_raw_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596 4.5M Nov 17 2023 5705STDY8058285_web_summary.html
-rw-r--r-- 1 addisony ds596 2.8M Nov 17 2023 cell_annotation.csv
-rw-r--r-- 1 addisony ds596 7.0K Nov 17 2023 E-MTAB-11115.idf.txt
-rw-r--r-- 1 addisony ds596  38K Nov 17 2023 E-MTAB-11115.sdrf.txt
```
Total size confirmed: du -sh mousescRNAseq/ → 2.8G
Files used going forward:

*_filtered_feature_bc_matrix.h5 — one per sample, 6 total
cell_annotation.csv — cell type labels for all 40,531 cells

Files not used:

*_molecule_info.h5 — raw molecule info, not needed
*_raw_feature_bc_matrix.h5 — unfiltered counts, not needed
*_web_summary.html — QC report for browser viewing only
*_metrics_summary.csv — summary statistics only
E-MTAB-11115.idf.txt / .sdrf.txt — ArrayExpress metadata only

### Dataset 2: Mouse Brain Visium (E-MTAB-11114)
```bash
wget -r --no-parent ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/114/E-MTAB-11114/Files/
```
**Output:** Not captured — command was run before logging began.

**Result confirmed by `ls -lh mouseVisium/`:**
```
total 5.5G
-rw-r--r-- 1 addisony ds596  9.2K Sep 27 2022 E-MTAB-11114.idf.txt
-rw-r--r-- 1 addisony ds596   31K Sep 27 2022 E-MTAB-11114.sdrf.txt
-rw-r--r-- 1 addisony ds596  754M Sep 27 2022 ST8059048_cloupe.cloupe
-rw-r--r-- 1 addisony ds596   18M Sep 27 2022 ST8059048_filtered_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596  971  Sep 27 2022 ST8059048_metrics_summary.csv
-rw-r--r-- 1 addisony ds596  218M Sep 27 2022 ST8059048_molecule_info.h5
-rw-r--r-- 1 addisony ds596   19M Sep 27 2022 ST8059048_raw_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596  9.9M Sep 27 2022 ST8059048_spatial.tar.gz
-rw-r--r-- 1 addisony ds596  9.0M Sep 27 2022 ST8059048_web_summary.html
-rw-r--r-- 1 addisony ds596  977M Sep 27 2022 ST8059049_cloupe.cloupe
-rw-r--r-- 1 addisony ds596   18M Sep 27 2022 ST8059049_filtered_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596  969  Sep 27 2022 ST8059049_metrics_summary.csv
-rw-r--r-- 1 addisony ds596  205M Sep 27 2022 ST8059049_molecule_info.h5
-rw-r--r-- 1 addisony ds596   19M Sep 27 2022 ST8059049_raw_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596   12M Sep 27 2022 ST8059049_spatial.tar.gz
-rw-r--r-- 1 addisony ds596  9.3M Sep 27 2022 ST8059049_web_summary.html
-rw-r--r-- 1 addisony ds596  973M Sep 27 2022 ST8059050_cloupe.cloupe
-rw-r--r-- 1 addisony ds596   19M Sep 27 2022 ST8059050_filtered_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596  969  Sep 27 2022 ST8059050_metrics_summary.csv
-rw-r--r-- 1 addisony ds596  237M Sep 27 2022 ST8059050_molecule_info.h5
-rw-r--r-- 1 addisony ds596   21M Sep 27 2022 ST8059050_raw_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596   12M Sep 27 2022 ST8059050_spatial.tar.gz
-rw-r--r-- 1 addisony ds596  9.4M Sep 27 2022 ST8059050_web_summary.html
-rw-r--r-- 1 addisony ds596  716M Sep 27 2022 ST8059051_cloupe.cloupe
-rw-r--r-- 1 addisony ds596   14M Sep 27 2022 ST8059051_filtered_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596  974  Sep 27 2022 ST8059051_metrics_summary.csv
-rw-r--r-- 1 addisony ds596  205M Sep 27 2022 ST8059051_molecule_info.h5
-rw-r--r-- 1 addisony ds596   18M Sep 27 2022 ST8059051_raw_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596  9.8M Sep 27 2022 ST8059051_spatial.tar.gz
-rw-r--r-- 1 addisony ds596  8.7M Sep 27 2022 ST8059051_web_summary.html
-rw-r--r-- 1 addisony ds596  760M Sep 27 2022 ST8059052_cloupe.cloupe
-rw-r--r-- 1 addisony ds596   17M Sep 27 2022 ST8059052_filtered_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596  971  Sep 27 2022 ST8059052_metrics_summary.csv
-rw-r--r-- 1 addisony ds596  234M Sep 27 2022 ST8059052_molecule_info.h5
-rw-r--r-- 1 addisony ds596   19M Sep 27 2022 ST8059052_raw_feature_bc_matrix.h5
-rw-r--r-- 1 addisony ds596  9.7M Sep 27 2022 ST8059052_spatial.tar.gz
-rw-r--r-- 1 addisony ds596  8.6M Sep 27 2022 ST8059052_web_summary.html
```

Total size confirmed: du -sh mouseVisium/ → 5.5G
Files used going forward:

*_filtered_feature_bc_matrix.h5 — one per section, 5 total
*_spatial/spatial/ — extracted spatial folders, one per section, 5 total

Files not used:

*.cloupe — for 10x Loupe browser only, ~750MB each, ignore
*_molecule_info.h5 — raw molecule info, not needed
*_raw_feature_bc_matrix.h5 — unfiltered counts, not needed
*_web_summary.html — QC report for browser viewing only
*_metrics_summary.csv — summary statistics only
*.tar.gz — already extracted, originals can stay
E-MTAB-11114.idf.txt / .sdrf.txt — ArrayExpress metadata only

### 1. Extract Visium Spatial Folders

```bash
cd mouseVisium/
for f in *.tar.gz; do
    sample="${f%%_spatial.tar.gz}"
    mkdir -p "${sample}_spatial"
    tar -xzf "$f" -C "${sample}_spatial"
done
```
**Output:** No errors. Created 5 spatial folders:
- `ST8059048_spatial/spatial/`
- `ST8059049_spatial/spatial/`
- `ST8059050_spatial/spatial/`
- `S8059051_spatial/spatial/`
- `ST8059052_spatial/spatial/`

Each containing:
```
aligned_fiducials.jpg
detected_tissue_image.jpg
scalefactors_json.json
tissue_hires_image.png
tissue_lowres_image.png
tissue_positions_list.csv
```

### 2. Pull Singularity Container
```bash
singularity pull cell2location.sif docker://quay.io/vitkl/cell2location:latest
```
Output: Successfully downloaded and converted ~4.8GB of OCI blobs to SIF format. Container saved as `cell2location.sif`

### 3. Move Container to Data Root
```bash
mv mouseVisium/cell2location.sif .

```
**Output:** No errors.
```
data/
├── cell2location.sif
├── mousescRNAseq/
└── mouseVisium/
```

### 4. Verify Container and Imports
```bash
singularity exec cell2location.sif python -c "
import cell2location
import scanpy
import torch
print('cell2location version:', cell2location.__version__)
print('scanpy version:', scanpy.__version__)
print('GPU available:', torch.cuda.is_available())
"
```
**Output:**
```
AttributeError: module 'cell2location' has no attribute '__version__'
```
Note: This is harmless — cell2location and scanpy both imported successfully. The error is only because this version does not expose a __version__ attribute. torch also imported fine.

### 5. Inspect Annotation File
bash
head -5 mousescRNAseq/cell_annotation.csv
wc -l mousescRNAseq/cell_annotation.csv
head -1 mousescRNAseq/cell_annotation.csv
```
**Output:**
```
Cell ID,sample,annotation_1,annotation_1_print
5705STDY8058283_AAACCCAAGCCTATTG-1,5705STDY8058283,Ext_L23,22_Ext_L23
5705STDY8058283_AAACCCAAGGTCATAA-1,5705STDY8058283,Oligo_2,59_Oligo_2
5705STDY8058283_AAACCCACAACCCTCT-1,5705STDY8058283,Oligo_2,59_Oligo_2
5705STDY8058283_AAACCCAGTGCTATTG-1,5705STDY8058283,Astro_THAL_med,10_Astro_THAL_med

40532 mousescRNAseq/cell_annotation.csv  (40531 cells + 1 header)

Cell ID,sample,annotation_1,annotation_1_print
```

**Key facts confirmed:**
- 40,531 annotated cells — matches paper exactly
- Cell ID format is `{sampleID}_{barcode}` e.g. `5705STDY8058283_AAACCCAAGCCTATTG-1`
- `annotation_1` is the clean cell type label e.g. `Ext_L23`
- `annotation_1_print` has a number prefix for ordering e.g. `22_Ext_L23`

---

### 6. Confirmed Data Structure
```
data/
├── cell2location.sif
├── mousescRNAseq/
│   ├── 5705STDY8058280_filtered_feature_bc_matrix.h5
│   ├── 5705STDY8058280_molecule_info.h5
│   ├── 5705STDY8058280_raw_feature_bc_matrix.h5
│   ├── 5705STDY8058281_filtered_feature_bc_matrix.h5
│   ├── 5705STDY8058281_molecule_info.h5
│   ├── 5705STDY8058281_raw_feature_bc_matrix.h5
│   ├── 5705STDY8058282_filtered_feature_bc_matrix.h5
│   ├── ... (same pattern for 8058283, 8058284, 8058285)
│   ├── cell_annotation.csv
│   ├── E-MTAB-11115.idf.txt
│   └── E-MTAB-11115.sdrf.txt
└── mouseVisium/
    ├── ST8059048_filtered_feature_bc_matrix.h5
    ├── ST8059048_spatial/spatial/
    ├── ST8059049_filtered_feature_bc_matrix.h5
    ├── ST8059049_spatial/spatial/
    ├── ST8059050_filtered_feature_bc_matrix.h5
    ├── ST8059050_spatial/spatial/
    ├── ST8059051_filtered_feature_bc_matrix.h5
    ├── ST8059051_spatial/spatial/
    ├── ST8059052_filtered_feature_bc_matrix.h5
    ├── ST8059052_spatial/spatial/
    ├── E-MTAB-11114.idf.txt
    └── E-MTAB-11114.sdrf.txt
```

---

### Current Working Directory
```
/projectnb/ds596/projects/Team 7/Cell2Location-Project-Research-Reproduction-Extention

### create conda environment
```bash
conda env create -f project_environment.yml
```

### Update conda environment
```bash
conda env update --file project_environment.yml
```