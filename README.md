# cattlePanSVimp

Workflows for imputing pangenome-derived **structural variants** into cattle
genotype data, and for running GWAS on the imputed panels.

Structural variants are largely invisible to the SNP arrays used in routine dairy
genotyping. These scripts build combined SV+SNP reference panels from a
Minigraph-Cactus pangenome, measure how accurately SVs can be imputed from four
marker densities (LD chip, HD chip, WGS SNPs, RNA-seq SNPs), and then apply the
resulting panels to association analysis.

The SV panel itself is built in a companion repository,
[cattleHolPanSV](https://github.com/xyxss/cattleHolPanSV).

> **Scope.** This repository is code only. It was written for a Slurm HPC cluster,
> and reproducing it end to end needs the input datasets listed under
> [Data availability](#data-availability) plus substantial compute. Individual
> stages are readable and reusable on their own.

---

## Repository layout

```
scripts/
  1.imputation_test/     Benchmark: how accurately can SVs be imputed?
  2.run_imputation/      Production: impute SVs into large cohorts
  3.GWAS_afImp/          GWAS on the imputed genotypes
  4.GWAS_valid_cdcb173WGS/  Validation GWAS against a 173-animal WGS set
  sample.group.csv       Public accessions and breed labels for panel animals
docs/
  steps/                 Step-by-step walkthrough of individual scripts
config.sh.example        Site paths - copy to config.sh and edit
environment.yml          Conda environment
```

Each stage directory has a `00.codes/` subdirectory holding the AWK helpers and
inner worker scripts that the numbered top-level scripts submit as jobs.

---

## Setup

### 1. Conda environment

```bash
conda env create -f environment.yml
conda activate cattle-pan-sv-imp
```

### 2. Site configuration

All cluster-specific paths live in one file, which is git-ignored:

```bash
cp config.sh.example config.sh
$EDITOR config.sh
source config.sh
```

Every script asserts the variables it needs at the top and exits immediately
with the missing variable's name if you haven't sourced it. Nothing writes to a
default or guessed location.

### 3. Tools not available through conda

Install these yourself and point `SOFTWARE_DIR` at them. Versions are the ones
used for the published run:

| Tool | Version | Used for |
|---|---|---|
| Beagle | `beagle.06Aug24.a91.jar` | Phasing and imputation (primary) |
| Minimac4 | 4.1.6 | Imputation (comparison arm) |
| GCTA | 1.94.1 | Mixed-model GWAS, COJO |
| SLEMM | 0.90.1 | Mixed-model GWAS at scale |
| EMMAX | (release used for validation GWAS) | Validation GWAS |
| RTG Tools | — | `rtg vcfeval` concordance scoring |

### 4. Reference genome

The workflows use **ARS-UCD2.0**, as `${REF_DIR}/ARS_UCD_v2.0.fa`.

You also need, alongside it:

- `ARS_UCD_v2.0.fa.fai`
- `ARS_UCD_v2.0.sdf` — RTG-format reference, required by `rtg vcfeval`
- `ARS_UCD_v2.0.ref_repeat/rm.<CLASS>.bed` — RepeatMasker classes, for
  repeat-stratified accuracy
- `ARS_UCD_v2.0.ref_gff/gff.<FEATURE>.bed` — genic features, for
  annotation-stratified accuracy

> **A naming caveat that trips people up:** inside the Minigraph-Cactus seqfile
> the reference haplotype is *labelled* `bosTau9`, and that label propagates into
> graph paths and some VCF sample columns. It is only a label. The underlying
> assembly is ARS-UCD2.0 throughout. (By the usual UCSC convention `bosTau9`
> means ARS-UCD1.2, which is **not** what is used here.)

### Input files you have to supply

Nothing in this repository produces these. Each is a variable in
`config.sh.example`; the scripts abort naming the variable if it is unset, but
they cannot tell you the file is missing until the step that reads it runs.

| Variable | What it is | Where it comes from |
|---|---|---|
| `REF_FA` | ARS-UCD2.0 FASTA + `.fai` | [NCBI GCF_002263795.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002263795.3/) |
| `REF_RM_DIR` | `rm.<CLASS>.bed` per repeat class | RepeatMasker over `REF_FA` |
| `REF_CHR_MAP` | RefSeq accession → chromosome name | Built from the assembly report |
| `SAMPLE_GROUP_TABLE` | `<sample_id>\t<group>`, no header | Your cohort. `scripts/sample.group.csv` is the published one, in CSV form |
| `PANGENIE_VCF` | PanGenie graph genotypes, merged biallelic | Companion repo [cattleHolPanSV](https://github.com/xyxss/cattleHolPanSV), step 6 |
| `DEEPVARIANT_VCF` | DeepVariant SNV calls, same cohort | Companion repo, step 4.2 |
| `BOVINEHD_MANIFEST` | `BovineHD_B1.csv` | Illumina (registration required) |
| `PHENO_DIR` | Phenotypes and PC covariates | **CDCB, restricted access — not redistributable** |

The group labels in `SAMPLE_GROUP_TABLE` are matched by string in
`1.imputation_prepare_data.sh` (`Holstein`, `Holstein-X-Jersey`, `Jersey`). If
yours are spelled differently the cohorts come out empty rather than erroring.

`SOFTWARE_DIR` additionally needs Beagle, Minimac4, Eagle, GCTA, EMMAX, and
SLEMM, none of which install through `environment.yml` — see the comments in
`config.sh.example` for versions and sources.

---

## Running the workflows

### Stage 1 — Imputation benchmark

Establishes how well SV imputation works before trusting it on real cohorts.

```bash
cd scripts/1.imputation_test
bash 1.imputation_prepare_data.sh     # filter panel, define cohorts, split by chromosome
bash 2.split_data_Submit.sh           # cross-validation folds
```

Then run one or more marker-density scenarios. Each submits an array of Slurm
jobs that call `00.codes/imp_all.sh`:

```bash
bash 3.run_for_HDlevel_related.sh     # HD chip density
bash 3.run_for_LDlevel_related.sh     # LD chip density
bash 3.run_for_WGSlevel_related.sh    # WGS SNP density
bash 3.run_for_RNAlevel_related.sh    # RNA-seq SNP density
```

`imp_all.sh` splits training and validation samples, builds a Beagle reference,
masks genotypes at set missingness rates, imputes, and scores the result with
`rtg vcfeval` — broken out by genotype class, SV type, repeat class, and genic
feature. Results accumulate in `<comp>.chr<N>.check.table2.GT.txt`.

Optional QC:

```bash
bash 4.PCA.sh          # PCA and Fst, to check population structure
bash 5.cds.check.sh    # coding-sequence overlap checks
```

### Stage 2 — Production imputation

```bash
cd scripts/2.run_imputation
bash 00.codes/0.refImp.sh   # build and index the Beagle reference panel
bash run_imp.sh             # chunked imputation driver
bash sub_imp.sh             # Slurm array submission, with job dependencies
```

`run_imp.sh` is the driver. It splits the reference and target by chromosome,
splits samples into chunks, runs Beagle per chunk, then merges back across
chunks and chromosomes. That chunking is what makes cohorts of tens of thousands
of animals tractable.

`sub_imp.sh` is the submission layer: it generates `sbatch` array scripts and
wires up `--dependency=afterok` so merging only starts once imputation finishes.

### Stage 3 — GWAS on imputed genotypes

```bash
cd scripts/3.GWAS_afImp
bash 1.prepare_genotypes.sh    # imputed VCF -> PLINK, filtered by MAF/R2
bash 2.prepare_phenotypes.sh   # phenotypes and covariates
bash run_slemm_mlm.sh          # or run_gcta_gwas.sh / run_plink_gwas.sh
```

Four association engines are included so results can be compared across methods:
SLEMM, GCTA (MLMA and fastGWA), PLINK, and EMMAX.

### Stage 4 — Validation GWAS

```bash
cd scripts/4.GWAS_valid_cdcb173WGS
bash data_pre.sh
bash GWAS_emmax_cdcb.sh
```

Repeats the association analysis in a 173-animal sequenced set, to check that
signals found on imputed genotypes hold up against directly sequenced ones.

---

## Script index

### `1.imputation_test/`

| Script | Purpose |
|---|---|
| `1.imputation_prepare_data.sh` | Filter the pangenome VCF, split by variant class, define benchmark cohorts, export per chromosome |
| `2.split_data_Submit.sh` | Build cross-validation folds and subsets |
| `3.run_for_HDlevel_related.sh` | Submit benchmark runs at HD-chip marker density |
| `3.run_for_LDlevel_related.sh` | Submit benchmark runs at LD-chip marker density |
| `3.run_for_WGSlevel_related.sh` | Submit benchmark runs at WGS-SNP density |
| `3.run_for_RNAlevel_related.sh` | Submit benchmark runs using RNA-seq-derived SNPs |
| `4.PCA.sh` | PCA and Fst across breed groups |
| `5.cds.check.sh` | Coding-sequence overlap QC for imputed variants |
| `00.codes/imp_all.sh` | Core worker: split, mask, impute, score one cohort × chromosome |
| `00.codes/ld_sv.sh` | LD between SVs and nearby SNPs |
| `00.codes/*.awk` | Helpers: VCF ID assignment, genotype masking, ploidy conversion, LD summaries |

### `2.run_imputation/`

| Script | Purpose |
|---|---|
| `run_imp.sh` | Chunked imputation driver: per-chromosome and per-sample-chunk Beagle, then merge |
| `sub_imp.sh` | Slurm array submission with `afterok` dependencies between imputation and merge |
| `00.codes/0.refImp.sh` | Build and index the Beagle reference panel |
| `00.codes/varImp.sh` | Variant imputation driver |
| `00.codes/varImp_chunk.sh` | Per-chunk variant imputation |
| `00.codes/varImp_mchk.sh` | Merge across sample chunks |
| `00.codes/varImp_mchr.sh` | Merge across chromosomes |

### `3.GWAS_afImp/`

| Script | Purpose |
|---|---|
| `1.prepare_genotypes.sh` | Imputed VCF to PLINK, MAF and imputation-R² filtering |
| `2.prepare_phenotypes.sh` | Phenotype and covariate tables |
| `run_slemm_mlm.sh` | SLEMM mixed-model GWAS |
| `run_SLEMM_afImpData.sh` | SLEMM on the imputed panel |
| `run_SLEMM_CDCB173WGS_val.sh` | SLEMM on the 173-animal WGS validation set |
| `run_gcta_gwas.sh` | GCTA fastGWA |
| `run_gcta-mlma.sh` | GCTA MLMA |
| `run_plink_gwas.sh` | PLINK association (baseline) |
| `phcor.awk` | Phenotype correlation helper |
| `00.codes/varImp*.sh` | Imputation drivers reused by this stage |

### `4.GWAS_valid_cdcb173WGS/`

| Script | Purpose |
|---|---|
| `data_pre.sh` | Assemble the validation genotype and phenotype inputs |
| `GWAS_emmax_cdcb.sh` | EMMAX GWAS on the validation set |

---

## Reproducibility notes

**Random sampling is seeded.** Benchmark cohorts (`Holstein750`, `Holstein250`,
`HolsteinRelated250`, `Jersey250`, and the `MultiBreed750` set built from the
three 250-animal groups) are drawn with a seeded shuffle in
`1.imputation_prepare_data.sh`, controlled by `SAMPLING_SEED` (default `1`, the
value used for the published results). The genotype-masking AWK helpers use fixed
`srand()` seeds. Changing any of these changes which animals and which genotypes
enter the benchmark, and therefore the accuracy estimates.

One caveat: seeded `shuf` is reproducible for a given GNU coreutils version, not
guaranteed identical across major coreutils releases. Record your coreutils
version alongside the seed if exact reproduction matters.

**Version pinning.** `environment.yml` pins minimum versions. The externally
installed tools are listed with exact versions in [Setup](#3-tools-not-available-through-conda).

**Line endings.** All files are LF, enforced by `.gitattributes`. This matters:
CRLF endings make `bash` fail to parse `function name (){` and leave stray `\r`
characters in AWK fields and CSV columns.

You can check every script parses before submitting anything:

```bash
find scripts -name '*.sh' -exec bash -n {} \;
```

### Known issues

- `scripts/1.imputation_test/00.codes/imp_all.sh` calls `ref2minimac4` twice, in
  two consecutive identical `if` blocks. It is idempotent, so results are
  unaffected, but the second call is wasted work.
- `2.run_imputation/00.codes/` and `3.GWAS_afImp/00.codes/` hold near-duplicate
  copies of the `varImp*` drivers. `varImp.sh` is identical between them;
  `varImp_chunk.sh`, `varImp_mchk.sh`, and `varImp_mchr.sh` have diverged
  slightly. They were not merged here because the differences are real and
  unreviewed — check which copy you are editing.
- Every script now parses (`bash -n` clean). Three had unbalanced loops, since
  repaired — they had been run interactively, section by section, rather than
  executed top to bottom:
  - `run_gcta_gwas.sh` and `run_plink_gwas.sh` carried stray `done` lines. Each
    submission block opens one loop over traits and closed it twice; the
    structurally identical blocks elsewhere in the same files close it once.
  - `1.prepare_genotypes.sh` opened `seq 1 29 | while read chr; do` that was
    never closed. `$chr` is not referenced anywhere in the file, and the inner
    loop already iterates over the per-chromosome VCFs directly, so the outer
    loop was removed rather than closed. Closing it instead would have submitted
    every job 29 times.

---

## Data availability

This repository holds **code only**. It contains no genotypes, phenotypes, or
sequence data, and none should be added — see `CONTRIBUTING.md`.

- `scripts/sample.group.csv` lists public BioSample and BioProject accessions
  with breed labels for the panel animals. These are public NCBI identifiers.
- The imputation panel and other large files are distributed separately; see the
  data availability statement in the manuscript.
- CDCB-derived genotypes and phenotypes used in stages 3 and 4 are
  **restricted-access** and cannot be redistributed here. Access requires an
  agreement with the Council on Dairy Cattle Breeding.

---

## Citation

See `CITATION.cff`. If you use this code, please cite the accompanying
manuscript once it is published.

## License

MIT — see `LICENSE`.

## Contact

Maintainer contact is listed in `CITATION.cff`. For questions about the
workflows, please open an issue.
