# Step 2.0 — Prepare input VCFs and sample groups for SV imputation

This step corresponds to the shell script:

- `scripts/2.0.imputation_prepare_data.sh` (original HPC version)
- `scripts_sanitized/2.0.imputation_prepare_data.sh` (public/sanitized version)

## Goal

Create standardized, analysis-ready VCFs (split by variant class, filtered, optionally thinned), and generate sample lists for downstream:
- cross-validation / subsetting,
- SNP-density experiments (WGS/HD/LD/RNA),
- SV imputation runs.

## Inputs

### Reference files
- `ARS_UCD_v2.0.fa` (reference FASTA)
- repeat annotation file used in some filtering steps (e.g., `ARS_UCD_v2.0.ref_repeat`)

> The script expects these under `ref_path`. In the public version, replace with your own reference locations.

### Sample grouping table
- `holPub.samples.group` (tab-delimited)
  - Column 1: sample ID
  - Column 2: group label (e.g., `Holstein`, `Jersey`, `Holstein-X-Jersey`)

The script generates the following sample lists:
- `holPub.samples` (all)
- `Holstein` (all Holstein)
- `Holstein750` (random 750 Holstein)
- `HolsteinRelated250` (random 250 Holstein-X-Jersey)
- `Jersey250` (random 250 Jersey)

### VCFs (examples inferred from the script)
The script processes multiple VCF sources, notably:
- a PanGenie VCF (`holPub.pangenie.vcf.gz`)
- a DeepVariant/GATK-like SNP VCF (`holPub.deepv.filter.vcf.gz` in later blocks)
- a minigraph / graph-derived VCF (`*.anno_biallelic.vcf.gz`)

Exact filenames and paths are set inside the script; the public/sanitized version replaces lab-specific absolute paths with placeholders.

### Helper AWK scripts
This step relies on multiple helper scripts under `00.codes/`, for example:
- `pan.indel.awk`
- `pan.complex.awk`
- `pangenie_remove50.awk`
- `pangenie.addID.awk`

In your upload, these live in `1.shell/3.run_for_imputation/code/`. For a public repo, place them under `awk/` (or `code/awk/`) and update references accordingly.

## Outputs (typical)

Depending on which blocks you run, expected outputs include (names may vary):

- Variant-class VCFs:
  - `holPub.pangenie-indels*.vcf.gz`
  - `holPub.pangenie-complex*.vcf.gz`
  - `holPub.pangenie-snps*.vcf.gz`
  - `holPub.pangenie-sv*.vcf.gz` (SV subset / filtered SVs)
- Filtered VCFs (e.g., with allele count filters, length filters)
- Thinned VCFs for SNP-density experiments:
  - `*.thin1k.vcf.gz`, `*.thin5k.vcf.gz`, etc. (created by `vcftools --thin`)
- Tabix indexes (`*.tbi`) for each bgzipped VCF
- Sample list files described above

## How to run

### On a SLURM cluster
Example (adjust cpus/mem/time for your environment):

```bash
sbatch -c 8 --mem=32G -t 24:00:00 scripts_sanitized/2.0.imputation_prepare_data.sh
```

The script uses:
- `nthreads=$SLURM_CPUS_PER_TASK`

So if you run interactively, you should set `SLURM_CPUS_PER_TASK` manually or modify `nthreads`.

### Locally / interactive
You can run interactively if you have `bcftools`, `vcftools`, `bgzip/tabix` available:

```bash
export SLURM_CPUS_PER_TASK=8
bash scripts_sanitized/2.0.imputation_prepare_data.sh
```

## Notes and common gotchas

1. **Chromosome selection**
   Several `bcftools view -r` calls assume autosomes `1..29` for cattle.

2. **Background jobs**
   Some blocks run with trailing `&`. Make sure to `wait` or monitor jobs before downstream steps.

3. **Potential filename mismatch**
   One block indexes `holPub.pangenie-sv3.vcf.gz` after writing `holPub.pangenie-sv.vcf.gz`.
   If that is not a deliberate naming scheme, correct it to keep filenames consistent.

4. **Reproducibility of random sampling**
   The sample subsets use `shuf -n ...` without a fixed seed.
   If you need exact reproducibility, replace with a seeded approach (e.g., `shuf --random-source=...`).

## What this enables next
- Step `2.1.subdata.vx.sh` (creating SNP-density subsets and CV folds)
- Step `2.10/2.11/2.12/2.13` (HD/LD/WGS/RNA imputation)
- Step `4.1.pca.sh` (PCA/Fst/LD-based downstream analyses)
