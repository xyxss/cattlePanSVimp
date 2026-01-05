# Cattle pangenome SV imputation workflows

This repository contains the analysis scripts used for generating pangenome-based variant panels (SV and SNP subsets), running imputation experiments under multiple marker-density scenarios (LD / HD / WGS-SNP / RNA-SNP subsets), and downstream QC / population structure analyses (PCA, Fst, LD).

> **Note:** The original scripts were written for an HPC environment (Slurm) and include site-specific absolute paths.  

## What’s in this repo

- `scripts_from_user/` – original HPC scripts (as provided)
- `scripts_sanitized/` – the same scripts with paths replaced by placeholders
- `docs/manuscript_snapshot.docx` – manuscript snapshot for context
- `environment.yml` – Conda environment (bcftools/vcftools/plink2 + R)
- `CITATION.cff` – citation metadata (edit `repository-code` after you create the GitHub repo)

## Quick start

### 1) Create the environment
```bash
conda env create -f environment.yml
conda activate cattle-pan-sv-imp
```

### 2) Set project paths
Most scripts expect a working directory and a reference:
```bash
export PROJECT_ROOT=/path/to/your/project
export WORK_DIR=$PROJECT_ROOT/work
export REF_FA=$PROJECT_ROOT/ref/ARS_UCD_v2.0.fa
```

### 3) Run the core steps (high level)

1. **Prepare variant panels / subsets**  
   - `2.0.imputation_prepare_data.sh`
   - `2.1.subdata.vx.sh`

2. **Imputation experiments**  
   - LD-chip subset: `2.11.ld_imputation.sh`
   - HD-chip subset: `2.10.hd_imputation.sh`
   - WGS-SNP subset: `2.12.wgs_imputation.sh`
   - RNA-SNP subset: `2.13.rna.sh`

3. **Population structure / stratification checks**  
   - PCA + Fst: `4.1.pca.sh`
   - LD estimation workflows: `3.1.subdata.ld.sh`

Because these scripts are HPC-oriented, typical usage is via `sbatch ... --wrap="bash <script> ..."` as shown within each script.

## Data availability

The manuscript snapshot includes a data-availability statement pointing to a download portal for the imputation panel (see `docs/manuscript_snapshot.docx`, “Availability of data and material”).

This GitHub repository is intended to host **code only**. Large datasets (VCF/BAM/CRAM, panels, etc.) should be hosted externally (e.g., institutional portal, NCBI/ENA, Zenodo, figshare) and linked from the README.

## How to cite

Update `CITATION.cff` with the final repository URL and (once available) the manuscript DOI.

## Contact

- Maintainer: Liu Yang


## Step-by-step (script-by-script)

If you prefer to reproduce the workflow one script at a time, start here:

- `docs/steps/02_imputation_prepare_data.md`



holPub.samples.group == sample.group.csv

hol-pg2hic-2024-05-22_graph_genotyping.merge-biallelic.vcf.gz