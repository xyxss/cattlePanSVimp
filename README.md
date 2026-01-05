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


## Step-by-step (script-by-script)

If you prefer to reproduce the workflow one script at a time, start here:

- `docs/steps/02_imputation_prepare_data.md`


# Quickstart (example)

This is a **minimal example** showing how to run the sanitized pipeline on your cluster.

> The scripts in this repo are sanitized: you must replace placeholder paths and ensure required tools are available (module/conda/singularity).

## 0) Create a conda environment (example)

```bash
conda env create -f environment.yml
conda activate cattleHolPanSV
```

## 1) Configure paths

Most scripts expect variables like:

- `work_dir` (project working directory)
- `ref_fa` (reference FASTA)
- input VCF(s) (SV/SNP panels)
- sample groups (see `scripts/sample.group.csv`)

Edit the script you want to run and replace placeholders, for example:

- `/PATH/TO/PROJECT`
- `/PATH/TO/REF/bosTau9.fa`
- `/PATH/TO/VCF/input.vcf.gz`

## 2) Prepare imputation benchmark inputs

```bash
cd scripts/1.imputation_test
bash 1.imputation_prepare_data.sh
```

Expected outputs are written under your configured output directory and include:
- filtered / subset VCFs
- sample lists / folds
- tabix indexes

## 3) Run one benchmark scenario (HD as an example)

```bash
cd scripts/1.imputation_test
bash 3.run_for_HDlevel_related.sh
```

This submits multiple SLURM jobs (via `sbatch`) that call the core wrapper:
- `scripts/1.imputation_test/00.codes/imp_all.sh`

## 4) PCA / QC (optional)

```bash
cd scripts/1.imputation_test
bash 4.PCA.sh
```

## 5) GWAS on imputed data (optional)

```bash
cd scripts/3.GWAS_afImp
bash 1.prepare_genotypes.sh
bash 2.prepare_phenotypes.sh
bash run_gcta_gwas.sh
```

## Notes

- Many scripts run multiple jobs in parallel; check your scheduler limits.
- For full reproducibility, pin tool versions (see `environment.yml`) and record seed values when sampling with `shuf`.


# Script index (sanitized)

This repository contains sanitized SLURM/bash pipelines used in the project.

## Layout
- `1.imputation_test/`: benchmark workflows (data prep, subset splits, HD/LD/WGS/RNA runs, PCA/QC)
- `2.run_imputation/`: production/batch imputation runners
- `3.GWAS_afImp/`: GWAS using imputed datasets (SLEMM/GCTA/PLINK)
- `4.GWAS_valid_cdcb173WGS/`: validation GWAS on CDCB 173 WGS

## Scripts

| Script | Purpose |
|---|---|
| `1.imputation_test/00.codes/imp_all.sh` | Core wrapper to run PanGenie-based genotyping/imputation on a given sample set and variant class. |
| `1.imputation_test/00.codes/ld_sv.sh` | Compute LD metrics for SV/SNP panels (used for QC/summary). |
| `1.imputation_test/1.imputation_prepare_data.sh` | Prepare VCF inputs, sample groups, and derived datasets for imputation benchmarks. |
| `1.imputation_test/2.13.rna.sh` | Imputation/genotyping evaluation using RNA-seq–derived variants or RNA-based sample sets. |
| `1.imputation_test/2.split_data_Submit.sh` | Split data into folds / subsets for benchmarking. |
| `1.imputation_test/3.run_for_HDlevel_related.sh` | Submit/coordinate imputation benchmark runs for HD-level chip-like density. |
| `1.imputation_test/3.run_for_LDlevel_related.sh` | Submit/coordinate imputation benchmark runs for LD-level chip-like density. |
| `1.imputation_test/3.run_for_WGSlevel_related.sh` | Submit/coordinate imputation benchmark runs for WGS-level density. |
| `1.imputation_test/4.PCA.sh` | Run PCA on selected genotype matrices / VCF-derived datasets. |
| `1.imputation_test/5.cds.check.sh` | Coding-sequence or CDS-related checks on imputed variants (QC). |
| `2.run_imputation/00.codes/0.refImp.sh` | Build reference indices/panels for imputation. |
| `2.run_imputation/00.codes/run.imp.sh` | (to be documented) |
| `2.run_imputation/00.codes/varImp.sh` | Variant imputation pipeline driver (chunked / multi-chrom / checks). |
| `2.run_imputation/00.codes/varImp_chunk.sh` | Variant imputation pipeline driver (chunked / multi-chrom / checks). |
| `2.run_imputation/00.codes/varImp_mchk.sh` | Variant imputation pipeline driver (chunked / multi-chrom / checks). |
| `2.run_imputation/00.codes/varImp_mchr.sh` | Variant imputation pipeline driver (chunked / multi-chrom / checks). |
| `2.run_imputation/run_imp.sh` | Batch submission wrappers for imputation jobs (Slurm). |
| `2.run_imputation/sub_imp.sh` | Batch submission wrappers for imputation jobs (Slurm). |
| `3.GWAS_afImp/00.codes/varImp.sh` | Variant imputation pipeline driver (chunked / multi-chrom / checks). |
| `3.GWAS_afImp/00.codes/varImp_chunk.sh` | Variant imputation pipeline driver (chunked / multi-chrom / checks). |
| `3.GWAS_afImp/00.codes/varImp_mchk.sh` | Variant imputation pipeline driver (chunked / multi-chrom / checks). |
| `3.GWAS_afImp/00.codes/varImp_mchr.sh` | Variant imputation pipeline driver (chunked / multi-chrom / checks). |
| `3.GWAS_afImp/1.prepare_genotypes.sh` | Convert imputed VCF/genotypes into GWAS-ready formats (plink/TPED, etc.). |
| `3.GWAS_afImp/2.prepare_phenotypes.sh` | Clean/format phenotype tables and covariates for GWAS. |
| `3.GWAS_afImp/run_SLEMM_CDCB173WGS_val.sh` | Run SLEMM mixed model GWAS / related analyses. |
| `3.GWAS_afImp/run_SLEMM_afImpData.sh` | Run SLEMM mixed model GWAS / related analyses. |
| `3.GWAS_afImp/run_gcta-mlma.sh` | Run GCTA GWAS/MLMA pipelines. |
| `3.GWAS_afImp/run_gcta_gwas.sh` | Run GCTA GWAS/MLMA pipelines. |
| `3.GWAS_afImp/run_plink_gwas.sh` | Run PLINK GWAS pipeline (baseline). |
| `3.GWAS_afImp/run_slemm_mlm.sh` | Run SLEMM mixed model GWAS / related analyses. |
| `4.GWAS_valid_cdcb173WGS/GWAS_emmax_cdcb.sh` | Run EMMAX GWAS (validation on CDCB 173 WGS). |
| `4.GWAS_valid_cdcb173WGS/data_pre.sh` | Prepare validation dataset inputs for GWAS. |



## Data availability

The manuscript snapshot includes a data-availability statement pointing to a download portal for the imputation panel (see `docs/manuscript_snapshot.docx`, “Availability of data and material”).

This GitHub repository is intended to host **code only**. Large datasets (VCF/BAM/CRAM, panels, etc.) should be hosted externally (e.g., institutional portal, NCBI/ENA, Zenodo, figshare) and linked from the README.

## How to cite

Update `CITATION.cff` with the final repository URL and (once available) the manuscript DOI.

## Contact

- Maintainer: Liu Yang
