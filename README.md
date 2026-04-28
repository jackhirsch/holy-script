# holy-script
An automated GWAS summary statistic munging script. Automatically cleans, detects columns, merges RSIDs, and prepares GWAS summary statistics for downstream analysis.

The final output is a tab-delimted file with the following columns:
RSID  CHR  POS  A1  A2  MAF  BETA  SE  Z  P  N

Script can handle gzipped files and files with any columns.
If it is impossible to generate one of the final columns, it will be assigned as NA.


To run, simply ./holy-script "sumstat_file" "genome build # (37 or 38)" "ancestry (EUR/AFR/AMR/EAS)" "N (N total or Neff works)"
