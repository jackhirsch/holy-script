<img 
  src="https://github.com/user-attachments/assets/82ffe31c-52d5-482a-b4fb-86d00cab7d9a"
  style="width: 25%; height: auto; float: left; margin-right: 15px;"
  alt="Picture2"
/>

# Holy-script for GWAS - GWAS Summary Statistics Harmonizer
An automated GWAS summary statistic munging script. Automatically cleans, detects columns, merges RSIDs, and prepares GWAS summary statistics for downstream analysis.

The final output is a tab-delimted file with the following columns:
RSID  CHR  POS  A1  A2  MAF  BETA  SE  Z  P  N

Script can handle gzipped files and files with any columns.

If it is impossible to generate one of the final columns, it will be assigned as NA.


To run, simply ./holy-script "sumstat_file" "genome build # (37 or 38)" "ancestry (EUR/AFR/AMR/EAS)" "N (N total or Neff works)"

## Table of Contents

- [Overview](#overview)
- [Output Format](#output-format)
- [Installation](#installation)
- [Usage](#usage)
- [Arguments](#arguments)
- [Pipeline Steps](#pipeline-steps)
- [Reference Files](#reference-files)
- [Statistical Notes](#statistical-notes)
- [Logging](#logging)
- [Examples](#examples)
- [FAQ](#faq)
- [License](#license)

---

## Overview

```
Input GWAS sumstats  ──►  holy-script  ──►  Harmonized TSV
(.txt/.tsv/.csv/.gz)                        (RSID CHR POS A1 A2 MAF BETA SE Z P N)
```

**What it does:**

1. Auto-detects file format and delimiter (tab / comma / space / gzip)
2. Fuzzy-matches column names for all common GWAS software output styles
3. Recovers missing RSIDs from CHR:POS:A2:A1 keys against a chromosome-split reference
4. Recovers missing CHR/POS from a RSID → position lookup table
5. Merges 1000 Genomes EUR allele frequencies
6. Aligns alleles (handles strand complements)
7. Corrects MAF > 0.5 (flips alleles, Z, beta)
8. Derives Z from BETA/SE, or from P-value, or approximates from MAF+N
9. Reconstructs BETA/SE when missing
10. Writes a clean, tab-separated output with `NA` for all missing values

---

## Output Format

Tab-separated, one row per SNP:

| Column | Type    | Description                          |
|--------|---------|--------------------------------------|
| RSID   | string  | rs identifier (e.g. rs12345)         |
| CHR    | string  | Chromosome (1–22, X, Y, MT)          |
| POS    | integer | Base-pair position (build 37 or 38)  |
| A1     | string  | Effect allele (always ≤ 0.5 MAF)     |
| A2     | string  | Other allele                         |
| MAF    | float   | Minor allele frequency (0–0.5)       |
| BETA   | float   | Effect size (log-OR for binary)      |
| SE     | float   | Standard error of BETA               |
| Z      | float   | Z-score                              |
| P      | float   | P-value (two-sided)                  |
| N      | float   | Sample size                          |

Missing values are written as `NA`.

---

## Installation

### Requirements

- R ≥ 4.1
- Packages: `data.table`, `optparse`
- [Git LFS](https://git-lfs.com) (for the bundled reference files)

```r
install.packages(c("data.table", "optparse"))
```

### Clone + setup

```bash
# 1. Install Git LFS (once per machine)
#    macOS:   brew install git-lfs
#    Ubuntu:  sudo apt install git-lfs
#    Conda:   conda install -c conda-forge git-lfs

or in cluster: module load git-lfs

# 2. Clone the repo
git clone https://github.com/jackhirsch/holy-script
cd holy-script

# 3. Pull the large reference files via LFS
bash setup_lfs.sh
```

After `setup_lfs.sh` completes, your folder will look like:

```
holy-script/
├── holy_script_gwas.R
├── README.md
├── LICENSE
├── .gitattributes
├── .gitignore
├── setup_lfs.sh
└── refs/
    ├── EA_freq.afreq              ← 1000G EUR MAF reference
    ├── build37_chr1.txt           ← RSID lookup, chr 1
    ├── build37_chr2.txt
    │   ...
    ├── build37_chr22.txt          ← RSID lookup, chr 22
    └── rsid_chrpos_build37.tsv    ← RSID → CHR/POS map
```

The reference files are tracked with **Git LFS** because they are too large for standard git. If you skip `setup_lfs.sh`, the files will be present as small LFS pointer stubs and the script will warn that they cannot be read.

---

## Usage

After cloning and running `setup_lfs.sh`, the reference paths default to `refs/` automatically — you don't need to pass them:

```bash
# Minimal — uses bundled refs/ by default
Rscript holy_script_gwas.R \
  --sumstat  my_gwas.tsv.gz \
  --build    37 \
  --ancestry EUR \
  --N        100000 \
  --output   my_gwas_harmonized.tsv
```

To override any reference path explicitly:

```bash
Rscript holy_script_gwas.R \
  --sumstat             my_gwas.tsv.gz \
  --build               37 \
  --ancestry            EUR \
  --N                   100000 \
  --maf_ref_path        /custom/path/EA_freq.afreq \
  --rsid_ref_path       /custom/path/build37_chr{CHR}.txt \
  --rsid_to_chrpos_path /custom/path/rsid_chrpos.tsv \
  --output              my_gwas_harmonized.tsv
```

---

## Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--sumstat` | ✅ | Input GWAS file (`.txt`, `.tsv`, `.csv`, `.gz`) |
| `--build` | ✅ | Genome build: `37` or `38` |
| `--ancestry` | ✅ | Ancestry (only `EUR` implemented) |
| `--N` | ⬜ | Override / supply sample size |
| `--maf_ref_path` | ⬜ | 1000G EUR `.afreq` file (default: `refs/EA_freq.afreq`) |
| `--rsid_ref_path` | ⬜ | Pattern for CHR:POS:A2:A1→RSID files, `{CHR}` = chromosome (default: `refs/build37_chr{CHR}.txt`) |
| `--rsid_to_chrpos_path` | ⬜ | RSID→CHR/POS mapping TSV (default: `refs/rsid_chrpos_build37.tsv`) |
| `--output` | ⬜ | Output path (default: `harmonized_gwas.tsv`) |

All reference arguments default to the bundled `refs/` folder. Override only if you want to use custom reference files.

---

## Pipeline Steps

### Step 1 — File Reading

- Auto-detects `.gz` (streaming) vs. plain text
- Sniffs delimiter from first line (tab, comma, space)
- Uses `data.table::fread()` for fast, memory-efficient loading
- Falls back to `sep = "auto"` on parse errors

### Step 2 — Column Detection (Fuzzy)

Detects the following columns by flexible name matching:

| Target | Accepted column names |
|--------|-----------------------|
| RSID | `rsid`, `rs_id`, `snp`, `snpid`, `markername`, `variantid`, `id` — or any column where ≥ 3 of the first 20 values start with `rs` |
| CHR | `chr`, `chrom`, `chromosome`, `#chrom` |
| POS | `pos`, `bp`, `position`, `basepair` |
| A1 (effect) | `a1`, `effect_allele`, `alt`, `allele1`, `tested_allele`, `ea` |
| A2 (other) | `a2`, `other_allele`, `ref`, `allele2`, `non_effect_allele`, `nea`, `oa` |
| BETA | `beta`, `logor`, `log_odds`, `effect`, `effect_size`, `b` |
| OR | `or`, `oddsratio`, `odds_ratio` → auto-converted to log(OR) |
| SE | any column containing `se` |
| Z | `z`, `zstat`, `z_score`, `zvalue`, `zstatistic` |
| P | `p`, `pval`, `pvalue`, `p_value`, `p.value` — or any column containing `p` |
| N | `n`, `nsamples`, `n_total`, `obs_ct`, `samplesize` |
| MAF | `maf`, `freq`, `af`, `eaf`, `a1_freq`, `effect_allele_frequency` |

Also handles combined `CHR:POS` columns (e.g. `1:12345`, `chr1_12345`).

### Step 3 — RSID Recovery

If RSID is missing but CHR + POS + A1 + A2 are present:

1. Constructs key: `CHR:POS:A2:A1`
2. Looks up per-chromosome reference files (`rsid_ref_path` pattern with `{CHR}`)
3. Uses `setkey` sorted join for efficiency

### Step 4 — CHR/POS Recovery

If CHR/POS are missing but RSID is present:

1. Looks up `rsid_to_chrpos_path`
2. Flexible column detection in the reference file itself

### Step 5 — 1000G MAF Merge

Reference file format (PLINK `.afreq`):

```
#CHROM  ID      REF  ALT  ALT_FREQS  OBS_CT
1       rs12345 A    G    0.1234     5008
```

Merged on RSID. ALT allele frequency is interpreted relative to the GWAS A1/A2 alleles.

### Step 6 — Allele Alignment

For each SNP with reference data:

```
A1 == ALT, A2 == REF  →  MAF = ALT_FREQ
A1 == REF, A2 == ALT  →  MAF = 1 - ALT_FREQ
Try complement strand  →  same rules on A/T C/G swaps
```

### Step 7 — MAF > 0.5 Correction

If MAF > 0.5 after alignment:

```
MAF  = 1 - MAF
A1  ↔  A2   (swap)
Z   = -Z    (flip sign)
BETA = -BETA (flip sign)
```

### Step 8 — Statistical Derivation

Priority order:

| Priority | Condition | Action |
|----------|-----------|--------|
| 1 | Z exists | Use directly |
| 2 | BETA + SE present | Z = BETA / SE |
| 3 | P exists | Z = Φ⁻¹(1 − P/2) × sign(BETA) |
| 4 | MAF + N + Z present | SE = 1/√(2·MAF·(1−MAF)·(N+Z²)); BETA = Z·SE ⚠️ approx |
| 5 | MAF + BETA + N only | SE = 1/√(2·N·MAF·(1−MAF)); Z = BETA/SE ⚠️ approx |
| — | Nothing available | NA |

All approximations are logged as `[WARNING]` with the formula applied.

---

## Reference Files

### MAF Reference (`--maf_ref_path`)

Generate from 1000 Genomes VCF using PLINK 2:

```bash
plink2 \
  --vcf ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
  --keep EUR_samples.txt \
  --freq \
  --out EA_freq
```

Output: `EA_freq.afreq`

### RSID Lookup (`--rsid_ref_path`)

One file per chromosome, tab-separated, no header:

```
CHR:POS:A2:A1<TAB>rs_identifier
```

Example filename pattern: `/ref/build37_chr{CHR}.txt`

Generate with (example awk from VCF):

```bash
for chr in {1..22}; do
  zcat vcf_chr${chr}.vcf.gz | awk -v chr=$chr '
    !/^#/ { print chr":"$2":"$4":"$5"\t"$3 }
  ' | LC_ALL=C sort -k1,1 > build37_chr${chr}.txt
done
```

### RSID → CHR/POS (`--rsid_to_chrpos_path`)

Tab-separated with header: `RSID CHR POS`

```
RSID    CHR  POS
rs12345 1    12345
```

---

## Statistical Notes

> **SE Approximation Warning**
>
> When BETA and SE are both missing, holy-script approximates SE using:
> ```
> SE = 1 / sqrt(2 · MAF · (1 − MAF) · (N + Z²))
> ```
> This assumes a linear genetic model and is an approximation. It is logged as `[WARNING]` and should be treated with caution in downstream analyses (e.g., LD Score regression).

> **P → Z Conversion**
>
> When Z is absent, it is derived from P as:
> ```
> Z = Φ⁻¹(P/2) × sign(BETA)
> ```
> If BETA is also absent, Z is assumed positive (sign = +1). This is noted in the log.

---

## Logging

Every run produces a `.log` file alongside the output (e.g., `my_gwas_harmonized.log`).

The log contains:
- Detected columns
- RSID recovery statistics
- MAF flip counts
- Allele swap counts
- Approximation warnings with formulas
- Final summary table

Example summary output:

```
[INFO] === HARMONIZATION SUMMARY ===
[INFO]   Input SNPs              : 8432190
[INFO]   Output SNPs             : 8432190
[INFO]   RSID recovered (pos key): 124301
[INFO]   CHR/POS recovered       : 0
[INFO]   Missing RSID (final)    : 3201
[INFO]   MAF flips               : 412890
[INFO]   Allele swaps            : 412890
[INFO]   Z computed from P       : 0
[INFO]   Beta/SE approximated    : 0
```

---

## Examples

### Standard run (bundled refs, default paths)

After `bash setup_lfs.sh`, just point at your sumstats:

```bash
Rscript holy_script_gwas.R \
  --sumstat  data/depression_gwas.tsv.gz \
  --build    37 \
  --ancestry EUR \
  --N        500000 \
  --output   out/depression_harmonized.tsv
```

### Override with custom reference files

```bash
Rscript holy_script_gwas.R \
  --sumstat             data/depression_gwas.tsv.gz \
  --build               37 \
  --ancestry            EUR \
  --N                   500000 \
  --maf_ref_path        /hpc/refs/EA_freq.afreq \
  --rsid_ref_path       /hpc/refs/build37_chr{CHR}.txt \
  --rsid_to_chrpos_path /hpc/refs/rsid_chrpos_build37.tsv \
  --output              out/depression_harmonized.tsv
```

### OR-based input (no BETA column)

holy-script auto-detects an `OR` column and converts to `log(OR)`:

```bash
Rscript holy_script_gwas.R \
  --sumstat or_based_gwas.tsv \
  --build 37 --ancestry EUR \
  --output harmonized.tsv
```

### Gzipped input, P-value only (no Z/BETA/SE)

Z will be derived from P; BETA/SE approximated from MAF+N:

```bash
Rscript holy_script_gwas.R \
  --sumstat pval_only.tsv.gz \
  --build 37 --ancestry EUR \
  --N 200000 \
  --maf_ref_path ref/EA_freq.afreq \
  --output harmonized.tsv
```

---

## FAQ

**Q: My file has no header row — will this work?**
A: `data.table::fread` with `header=TRUE` will treat the first row as a header. If your file truly has no header, add one before running. If an RSID column cannot be detected from the header, the script will scan column values for `rs`-prefixed entries.

**Q: What if I have both CHR/POS and RSID?**
A: RSID is used as the primary merge key. CHR/POS are preserved for output.

**Q: Does this work for binary traits (case/control)?**
A: Yes. Provide log(OR) as BETA or an OR column. The script converts OR → log(OR) automatically.

**Q: Can I run this on a SLURM cluster?**
A: Yes. The script has no interactive elements. Wrap it in a SLURM batch script:
```bash
#!/bin/bash
#SBATCH --mem=32G --cpus-per-task=4
module load R/4.3.0
Rscript holy_script_gwas.R --sumstat $SUMSTAT ...
```

**Q: What build is supported?**
A: Build 37 (hg19) and Build 38 (hg38). Provide matching reference files for your build.

**Q: Why is ancestry only EUR?**
A: 1000G MAF reference files for other populations (AFR, EAS, SAS, AMR) are not yet integrated but the architecture is designed to support them. PRs welcome.

---

## Contributing

Issues and pull requests are welcome. Please open an issue before starting large changes. This is a work in progress!

---

## License

MIT License. See [LICENSE](LICENSE).

---

## Citation (optional)

If you use holy-script in published research, feel free to cite this repository and acknowledge the 1000 Genomes Project for reference allele frequencies.
