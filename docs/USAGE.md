# Usage notes (mapping scripts to analyses)

## Variant panel preparation
- `2.0.imputation_prepare_data.sh`: filtering, grouping, SV-type splits, thinning, per-chromosome exports
- `2.1.subdata.vx.sh`: additional subsetting / comparisons

## Imputation scenarios
- `2.11.ld_imputation.sh`: SV + LD-SNP hybrid panels
- `2.10.hd_imputation.sh`: SV + HD-SNP hybrid panels; includes liftover notes (UMD3.1 ↔ ARS-UCD2.0)
- `2.12.wgs_imputation.sh`: SV + WGS-SNP hybrid panels
- `2.13.rna.sh`: SV + RNA-SNP hybrid panels

## Population structure & stratification
- `4.1.pca.sh`: PCA + Fst

## LD calculations
- `3.1.subdata.ld.sh`: LD computations and post-processing

## Portability
The scripts assume Slurm and reference helper AWK scripts (e.g. `haploid2diploid.awk`, `addID.awk`).
For a polished public release, consider adding those helper scripts under `awk/` and a `config/config.sh`.


## Script-by-script walkthrough

- **Step 2.0 — Prepare data for imputation**: `2.0.imputation_prepare_data.sh`
  - See: `docs/steps/02_imputation_prepare_data.md`
