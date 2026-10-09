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

A minimal repository can be organized as follows:

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

The `{CHR}` placeholder is replaced with the chromosome being processed. For example, chromosome 13 will use:

```text
/path/to/rsids_37/build37_allmultiallelic_chr13.txt
```

**Important:** supply a path that exists on your system. If the script has a machine-specific default path, explicitly setting `--rsid_ref_path` makes the command portable across computers and HPC systems.

For a custom frequency reference:

```bash
Rscript holy-script.R \
  --sumstat input_sumstats.tsv \
  --build 37 \
  --ancestry EUR \
  --maf_ref_path ref/EA_freq.afreq \
  --rsid_ref_path "/path/to/rsids_37/build37_allmultiallelic_chr{CHR}.txt" \
  --output holy_sumstats.tsv
```

## Reference files

`holy-script` can use two types of reference data:

1. Ancestry-specific allele-frequency references in PLINK 2 `.afreq` format.
2. Chromosome-specific mappings from `CHR:POS:REF:ALT` keys to rsIDs.

These files serve different purposes and are not interchangeable.

### 1. 1000 Genomes allele-frequency references

The frequency references provide allele frequencies that can be used to fill missing or invalid frequency information in the input summary statistics.

The expected default filenames are:

| Ancestry | Filename |
|---|---|
| European (EUR) | `ref/EA_freq.afreq` |
| African (AFR) | `ref/AFR_freq.afreq` |
| Admixed American (AMR) | `ref/AMR_freq.afreq` |

The EUR file is named `EA_freq.afreq` for compatibility with the script's default configuration.

#### Where to obtain the data

The 1000 Genomes Project Phase 3 reference panel is available through the International Genome Sample Resource (IGSR):

- 1000 Genomes Phase 3 data portal: https://www.internationalgenome.org/data-portal/data-collection/phase3
- 1000 Genomes FTP site: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
- PLINK 2 resources, including Phase 3 reference data: https://www.cog-genomics.org/plink/2.0/resources

For this workflow, use a reference panel with coordinates compatible with your GWAS data. The original Phase 3 release uses GRCh37 coordinates. If your summary statistics use GRCh38, use a matching GRCh38 reference instead.

#### Generating ancestry-specific `.afreq` files

The `.afreq` files can be generated from 1000 Genomes genotypes using PLINK 2.

1. Obtain the Phase 3 genotype files in PLINK 2 format, along with their sample metadata.
2. Identify the samples belonging to the desired superpopulation using the 1000 Genomes sample panel.
3. Create a sample-keep file for each ancestry.
4. Calculate allele frequencies for each ancestry with PLINK 2.

For example, once the genotype files and EUR sample-keep file are available:

```bash
plink2 \
  --pfile all_phase3 \
  --keep EUR.keep \
  --freq \
  --out EA_freq
```

Repeat the frequency calculation with the corresponding AFR and AMR sample lists:

```bash
plink2 \
  --pfile all_phase3 \
  --keep AFR.keep \
  --freq \
  --out AFR_freq

plink2 \
  --pfile all_phase3 \
  --keep AMR.keep \
  --freq \
  --out AMR_freq
```

The `--pfile` argument must match the prefix and format of the genotype files you downloaded. For compressed PLINK 2 files, use the corresponding PLINK 2 input options and filenames. The keep files must contain sample IDs that match the genotype dataset's sample IDs.

The resulting frequency files should be placed in the `ref/` directory or passed explicitly through `--maf_ref_path`.

**Note:** the commands above assume that the genotype dataset and ancestry-specific sample lists have already been prepared. They do not automatically download the genotype data or identify ancestry samples.

#### Required `.afreq` format

The script expects the standard PLINK 2 allele-frequency format with these columns:

```text
#CHROM	ID	REF	ALT	ALT_FREQS	OBS_CT
```

Example:

```text
#CHROM	ID	REF	ALT	ALT_FREQS	OBS_CT
7	rs10	C	A	0.00302572	1322
7	rs10000	A	G	0.0257186	1322
7	rs10001	C	T	0.157337	1322
```

The file is tab-separated and includes a header.

| Column | Meaning |
|---|---|
| `#CHROM` | Chromosome |
| `ID` | Variant identifier, usually an rsID |
| `REF` | Reference allele |
| `ALT` | Alternate allele |
| `ALT_FREQS` | Frequency of the alternate allele |
| `OBS_CT` | Number of allele observations used to calculate the frequency |

The values in `ALT_FREQS` are frequencies of the reference file's `ALT` allele. They are not automatically the frequency of the GWAS effect allele. Alleles must be aligned correctly before frequencies are compared or used.

Use the ancestry reference appropriate to the analysis. The script can skip loading this reference when the input already contains valid frequency information.

### 2. Chromosome-position-allele to rsID references

These files allow the script to recover rsIDs when they are missing from the GWAS summary statistics.

The file naming convention is:

```text
build37_allmultiallelic_chr1.txt
build37_allmultiallelic_chr2.txt
...
build37_allmultiallelic_chr22.txt
build37_allmultiallelic_chrX.txt
```

The chromosome placeholder in `--rsid_ref_path` must match this naming pattern.

#### Required file format

Each file must have **two columns, no header**, separated by whitespace:

1. A variant key in `CHR:POS:REF:ALT` format.
2. The corresponding rsID.

Example:

```text
13:19020002:A:G rs1870154548
13:19020002:A:AATT rs1870154627
13:19020006:C:G rs1593354050
13:19020007:A:C rs1870154750
13:19020007:A:T rs1870154750
13:19020010:T:A rs1870154811
13:19020012:C:T rs1396158199
13:19020013:C:A rs181615907
13:19020013:C:T rs181615907
13:19020014:TC:T rs1466826665
```

The columns are:

| Column | Meaning | Example |
|---|---|---|
| 1 | Chromosome, position, reference allele, and alternate allele joined with colons | `13:19020002:A:G` |
| 2 | rsID | `rs1870154548` |

Requirements:

- Do not add a header row.
- Use whitespace (tabs or spaces) between the two columns.
- Use the same chromosome naming convention as the script's parsed coordinates; the examples use `13`, not `chr13`.
- The first column must include both alleles, not just chromosome and position.
- Keep the chromosome in the filename consistent with the chromosome encoded in the first column.
- Use a reference file compatible with the genome build specified by `--build`.

**Why are both alleles required?**

Multiple variants can occur at the same genomic position. For example, a position can have multiple alternate alleles. Using `CHR:POS:REF:ALT` distinguishes these variants more precisely than chromosome and position alone.

The reference must use the same allele representation and genome build as the input data for the key to match.

#### Where to obtain or create these files

These files are lookup tables, not the standard 1000 Genomes `.afreq` files. They can be created from a suitable variant database or annotation source containing genomic coordinates, reference and alternate alleles, and rsIDs.

For GRCh37, NCBI's dbSNP resources and other variant annotation resources can be used to obtain variant identifiers and coordinates:

- NCBI dbSNP: https://www.ncbi.nlm.nih.gov/snp/
- NCBI Datasets: https://www.ncbi.nlm.nih.gov/datasets/

If starting with a downloaded variant annotation file, extract the chromosome, position, reference allele, alternate allele, and rsID fields, then format each record as:

```text
CHR:POS:REF:ALT rsID
```

Split the resulting table into chromosome-specific files. Ensure that the coordinates and alleles correspond to the intended genome assembly. A GRCh37 reference must not be treated as a GRCh38 reference simply by changing its filename.

For multiallelic sites, retain the appropriate allele-specific records. The same rsID may appear on multiple rows when multiple alternate alleles belong to that variant record or identifier.

### 3. Optional rsID-to-coordinate reference

The optional `--rsid_to_chrpos_path` argument accepts a mapping file that allows the script to recover genomic coordinates from rsIDs.

This file differs from the chromosome-specific reference files above.

#### Required format

A tab-separated file with a header containing columns for rsID, chromosome, and position:

```text
RSID	CHR	POS
rs123456789	1	10002
rs987654321	2	20003
```

Column names are matched flexibly, but the file must contain recognizable rsID, chromosome, and position columns.

Use this option when your input has rsIDs but lacks chromosome and position columns and you have an appropriate mapping reference available.

Example:

```bash
Rscript holy-script.R \
  --sumstat input_sumstats.tsv \
  --build 37 \
  --ancestry EUR \
  --rsid_to_chrpos_path ref/rsid_to_chrpos.tsv \
  --output holy_sumstats.tsv
```

## Input summary statistics

The script is designed to recognize common GWAS column names, including common representations of:

- Chromosome and base-pair position
- Variant ID or rsID
- Reference and alternate alleles
- Effect allele and other allele
- Effect estimate, standard error, and p-value
- Z statistic
- Sample size
- Allele frequency

For example, an input may contain columns such as:

```text
CHROM POS ID REF ALT A1 A1_FREQ TEST OBS_CT OR LOG(OR)_SE L95 U95 Z_STAT P ERRCODE
```

The exact columns required depend on the information present in the input and what can be recovered from the available references.

Before running the script, verify that:

- The genome build is known.
- Allele columns are interpreted correctly.
- Chromosome and position values are consistent with the specified build.
- Effect estimates and standard errors use compatible representations.
- The selected ancestry is appropriate for the frequency reference.

## Output

The script writes the holy summary-statistic file specified by `--output`.

It also writes a log file using the output filename with `.log` as the extension. For example:

```text
holy_sumstats.tsv
holy_sumstats.log
```

The log reports processing information and counters, including rsID recovery and output counts. Review the log to identify missing identifiers, unmatched variants, or other warnings.

## Troubleshooting

### No rsIDs were recovered

Check that:

- `--rsid_ref_path` points to the correct directory and filename pattern.
- The reference files have no header row.
- The first column uses `CHR:POS:REF:ALT` format.
- Chromosome labels and genome builds match the input.
- Reference and alternate alleles are represented consistently.

### The frequency reference is not loaded

If the input contains valid allele-frequency values, the script may skip loading the frequency reference. This is expected behavior.

If frequencies are missing or invalid, check that the reference file exists, has the expected PLINK 2 `.afreq` columns, and corresponds to the selected ancestry.

### Variants fail to match despite having the same position

A match can fail when the reference and input use different genome builds, allele representations, or reference/alternate allele orientations. Matching by position alone is not sufficient when multiple alleles occur at the same site.

### GRCh37 and GRCh38 references are mixed

Make sure the input summary statistics, chromosome-position-to-rsID lookup, and optional rsID-to-coordinate lookup use compatible assemblies. Do not use a GRCh37 lookup for GRCh38 data.

## Reproducibility and data management

Reference files can be large and may have their own licensing, attribution, or distribution requirements. Download or generate them separately rather than committing large reference datasets to this repository.

Record the source, genome build, ancestry, and generation method for each reference file. These details are important for reproducing harmonization results.

## References

- International Genome Sample Resource (IGSR), 1000 Genomes Project Phase 3: https://www.internationalgenome.org/data-portal/data-collection/phase3
- IGSR FTP site: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/
- PLINK 2 resources: https://www.cog-genomics.org/plink/2.0/resources
- NCBI dbSNP: https://www.ncbi.nlm.nih.gov/snp/
- NCBI Datasets: https://www.ncbi.nlm.nih.gov/datasets/
