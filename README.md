# holy-script

A command-line R script for harmonizing GWAS summary statistics into a consistent (holy) format for downstream genetic analyses.
Author: Jack Hirsch

`holy-script.R` reads GWAS summary statistics, identifies common column names, standardizes variant information, recovers missing rsIDs or genomic coordinates when reference files are provided, and uses ancestry-specific allele frequencies to fill missing frequency information when needed.

## Features

- Supports common GWAS summary-statistic formats and column names.
- Supports GRCh37 and GRCh38 coordinate builds.
- Supports EUR, AFR, and AMR ancestry-specific allele-frequency references.
- Recovers missing rsIDs using chromosome, position, reference allele, and alternate allele.
- Optionally recovers chromosome and position from rsIDs.
- Uses valid allele-frequency values already present in the input and only loads the ancestry frequency reference when needed. Ex. If standardizing an N-weighted meta-analysis and Beta and SE are not available, it can be stimated
- Z scores and Beta can be estimated from OR
- Writes holy summary statistics and a log file.

## Repository contents


```text
holy-script/
├── holy-script.R
├── README.md
└── ref/
    ├── EA_freq.afreq
    ├── AFR_freq.afreq
    └── AMR_freq.afreq
```

The `ref/` directory is optional if the input summary statistics already contain valid allele-frequency values. Reference files are not included in the repository and should be obtained or generated separately.

The chromosome-position-to-RSID reference files can be stored in a separate directory, particularly if they are large.

## Requirements

- R
- R packages:
  - `data.table`
  - `optparse`

Install the required packages in R:

```r
install.packages(c("data.table", "optparse"))
```

## Usage

Run the script from a terminal:

```bash
Rscript holy-script.R \
  --sumstat input_sumstats.tsv \
  --build 37 \
  --ancestry EUR \
  --output holy_sumstats.tsv
```

### Command-line arguments

| Argument | Description |
|---|---|
| `--sumstat` | Path to the input GWAS summary-statistic file. Required. |
| `--build` | Genome assembly: `37` or `38`. Defaults to `37`. |
| `--ancestry` | Ancestry for the frequency reference: `EUR`, `AFR`, or `AMR`. Defaults to `EUR`. |
| `--N` | Optional sample size to use when it cannot be obtained from the input. |
| `--maf_ref_path` | Optional path to an ancestry-specific PLINK 2 `.afreq` reference file. |
| `--rsid_ref_path` | Pattern for chromosome-specific `CHR:POS:REF:ALT` to rsID files. Use `{CHR}` as the chromosome placeholder. |
| `--rsid_to_chrpos_path` | Optional file mapping rsIDs to chromosome and position. |
| `--output` | Output filename. Defaults to `holy_gwas.tsv`. |

### Specifying reference files

For example, if your chromosome-specific rsID files are stored in `/path/to/rsids_37/`, use:

```bash
Rscript holy-script.R \
  --sumstat input_sumstats.tsv \
  --build 37 \
  --ancestry EUR \
  --rsid_ref_path "/path/to/rsids_37/build37_allmultiallelic_chr{CHR}.txt" \
  --output holy_sumstats.tsv
```

#### Where to obtain the data

The 1000 Genomes Project Phase 3 reference panel is available through the International Genome Sample Resource (IGSR):

- 1000 Genomes Phase 3 data portal: https://www.internationalgenome.org/data-portal/data-collection/phase3
- 1000 Genomes FTP site: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
- PLINK 2 resources, including Phase 3 reference data: https://www.cog-genomics.org/plink/2.0/resources
