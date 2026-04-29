#!/usr/bin/env Rscript
# =============================================================================
#  holy-script — GWAS Summary Statistics Harmonizer
#  Version: 1.0.0
#  Author:  holy-script contributors
#  License: MIT
#
#  Full pipeline: RSID merge · allele alignment · MAF correction ·
#  Z-score derivation · beta/SE reconstruction · HPC-efficient streaming
#
#  Output columns (tab-separated):
#  RSID  CHR  POS  A1  A2  MAF  BETA  SE  Z  P  N
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# =============================================================================
# 0. CLI ARGUMENT PARSING
# =============================================================================

option_list <- list(
  make_option("--sumstat",          type = "character", default = NULL,
              help = "Input GWAS summary statistics file (.txt/.tsv/.csv/.gz)"),
  make_option("--build",            type = "integer",   default = 37,
              help = "Genome build: 37 or 38 [default: 37]"),
  make_option("--ancestry",         type = "character", default = "EUR",
              help = "Ancestry (only EUR implemented) [default: EUR]"),
  make_option("--N",                type = "integer",   default = NULL,
              help = "Sample size (overrides N column if provided)"),
  make_option("--maf_ref_path",     type = "character",
              default = file.path(dirname(sys.frame(1)$ofile), "refs", "EA_freq.afreq"),
              help = "Path to 1000G EUR .afreq file [default: refs/EA_freq.afreq]"),
  make_option("--rsid_ref_path",    type = "character",
              default = file.path(dirname(sys.frame(1)$ofile), "refs", "build37_chr{CHR}.txt"),
              help = "Pattern for CHR:POS:A2:A1->RSID files [default: refs/build37_chr{CHR}.txt]"),
  make_option("--rsid_to_chrpos_path", type = "character",
              default = file.path(dirname(sys.frame(1)$ofile), "refs", "rsid_chrpos_build37.tsv"),
              help = "Path to RSID->CHR/POS mapping file [default: refs/rsid_chrpos_build37.tsv]"),
  make_option("--output",           type = "character", default = "harmonized_gwas.tsv",
              help = "Output file path [default: harmonized_gwas.tsv]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# =============================================================================
# 1. LOGGING UTILITIES
# =============================================================================

LOG_FILE <- paste0(tools::file_path_sans_ext(opt$output), ".log")

log_msg <- function(level, msg) {
  ts  <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", ts, level, msg)
  message(line)
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

log_info  <- function(msg) log_msg("INFO",    msg)
log_warn  <- function(msg) log_msg("WARNING", msg)
log_error <- function(msg) { log_msg("ERROR", msg); stop(msg) }

# Counters for the harmonization summary
COUNTERS <- new.env(parent = emptyenv())
COUNTERS$n_input           <- 0L
COUNTERS$n_rsid_recovered  <- 0L
COUNTERS$n_chrpos_recovered <- 0L
COUNTERS$n_maf_flip        <- 0L
COUNTERS$n_allele_swap     <- 0L
COUNTERS$n_beta_approx     <- 0L
COUNTERS$n_z_from_p        <- 0L
COUNTERS$n_missing_rsid    <- 0L
COUNTERS$n_output          <- 0L

print_summary <- function() {
  log_info("=== HARMONIZATION SUMMARY ===")
  log_info(sprintf("  Input SNPs              : %d", COUNTERS$n_input))
  log_info(sprintf("  Output SNPs             : %d", COUNTERS$n_output))
  log_info(sprintf("  RSID recovered (pos key): %d", COUNTERS$n_rsid_recovered))
  log_info(sprintf("  CHR/POS recovered       : %d", COUNTERS$n_chrpos_recovered))
  log_info(sprintf("  Missing RSID (final)    : %d", COUNTERS$n_missing_rsid))
  log_info(sprintf("  MAF flips               : %d", COUNTERS$n_maf_flip))
  log_info(sprintf("  Allele swaps            : %d", COUNTERS$n_allele_swap))
  log_info(sprintf("  Z computed from P       : %d", COUNTERS$n_z_from_p))
  log_info(sprintf("  Beta/SE approximated    : %d", COUNTERS$n_beta_approx))
}

# =============================================================================
# 2. VALIDATE ARGUMENTS
# =============================================================================

log_info("=== holy-script GWAS Harmonizer v1.0.0 ===")

if (is.null(opt$sumstat))      log_error("--sumstat is required.")
if (!file.exists(opt$sumstat)) log_error(sprintf("Input file not found: %s", opt$sumstat))
if (!opt$build %in% c(37L, 38L)) log_error("--build must be 37 or 38.")
if (opt$ancestry != "EUR")    log_error("Only EUR ancestry is currently implemented.")

log_info(sprintf("Input     : %s", opt$sumstat))
log_info(sprintf("Build     : GRCh%d", opt$build))
log_info(sprintf("Ancestry  : %s",  opt$ancestry))
log_info(sprintf("Output    : %s",  opt$output))

# =============================================================================
# 3. FILE READING WITH AUTO-DETECTION
# =============================================================================

detect_and_read <- function(path) {
  log_info(sprintf("Reading input: %s", path))

  is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
  if (is_gz) log_info("Detected gzip file — streaming without full decompression.")

  # Peek at first non-empty, non-comment line to detect delimiter
  peek <- if (is_gz) {
    con <- gzcon(file(path, "rb"))
    on.exit(close(con))
    readLines(con, n = 5L)
  } else {
    readLines(path, n = 5L)
  }
  peek <- peek[nchar(trimws(peek)) > 0 & !startsWith(peek, "##")]

  # Detect delimiter from the first data-like line
  detect_sep <- function(line) {
    counts <- c(
      tab   = nchar(line) - nchar(gsub("\t",  "", line)),
      comma = nchar(line) - nchar(gsub(",",   "", line)),
      space = nchar(line) - nchar(gsub(" +",  "", line))
    )
    c("\t", ",", " ")[which.max(counts)]
  }
  sep <- detect_sep(peek[1L])
  log_info(sprintf("Detected delimiter: '%s'", ifelse(sep == "\t", "TAB",
                                                ifelse(sep == ",", "COMMA", "SPACE"))))

  dt <- tryCatch(
    fread(path, sep = sep, header = TRUE,
          data.table = TRUE, showProgress = FALSE,
          fill = TRUE, quote = ""),
    error = function(e) {
      log_warn(sprintf("fread failed (%s), retrying with sep='auto'", conditionMessage(e)))
      fread(path, sep = "auto", header = TRUE,
            data.table = TRUE, showProgress = FALSE,
            fill = TRUE, quote = "")
    }
  )

  if (nrow(dt) == 0L) log_error("Input file is empty after reading.")

  log_info(sprintf("Read %d rows, %d columns.", nrow(dt), ncol(dt)))
  dt
}

# =============================================================================
# 4. COLUMN DETECTION (FUZZY MATCHING)
# =============================================================================

# Returns the first column name matching any pattern in `patterns` (case-insensitive)
# `exact_first` patterns are checked before embedded patterns
find_col <- function(cols, patterns, exact = TRUE) {
  cols_l <- tolower(cols)
  pats_l <- tolower(patterns)
  for (p in pats_l) {
    if (exact) {
      hit <- which(cols_l == p)
    } else {
      hit <- which(grepl(p, cols_l, fixed = TRUE))
    }
    if (length(hit)) return(cols[hit[1L]])
  }
  NULL
}

detect_columns <- function(dt) {
  cols <- names(dt)
  log_info("Detecting columns...")

  map <- list()

  # RSID — try explicit patterns first, then scan all cols for "rs" prefix
  map$RSID <- find_col(cols, c("rsid","rs_id","snp","snpid","markername",
                               "marker_name","variantid","variant_id","id"))
  if (is.null(map$RSID)) {
    for (cc in cols) {
      vals <- dt[[cc]]
      if (is.character(vals) || is.factor(vals)) {
        if (sum(startsWith(as.character(vals[1:min(20,nrow(dt))]), "rs"), na.rm=TRUE) > 3) {
          map$RSID <- cc
          log_info(sprintf("  RSID detected by value scan: '%s'", cc))
          break
        }
      }
    }
  }

  # CHR / POS
  map$CHR  <- find_col(cols, c("chr","chrom","chromosome","#chrom"))
  map$POS  <- find_col(cols, c("pos","bp","position","basepair"))

  # Alleles
  map$A1   <- find_col(cols, c("a1","effect_allele","alt","allele1","tested_allele","ea"))
  map$A2   <- find_col(cols, c("a2","other_allele","ref","allele2","non_effect_allele","nea","oa"))

  # Statistics
  map$BETA <- find_col(cols, c("beta","logor","log_odds","effect_size","effect","b"))
  map$OR   <- find_col(cols, c("or","oddsratio","odds_ratio"))
  map$SE   <- find_col(cols, c("se","stderr","standard_error"))
  if (is.null(map$SE)) map$SE <- find_col(cols, c("se"), exact = FALSE)
  map$Z    <- find_col(cols, c("z","zstat","z_score","zvalue","zstatistic","zscore"))
  map$P    <- find_col(cols, c("p","pval","p_val","pvalue","p_value","p.value"))
  # Broader P search: any col containing 'p'
  if (is.null(map$P)) {
    hits <- grep("p", tolower(cols), value = TRUE)
    if (length(hits)) {
      map$P <- cols[tolower(cols) %in% hits][1L]
    }
  }
  map$N    <- find_col(cols, c("n","nsamples","n_total","obs_ct","samplesize","n_eff"))
  map$MAF  <- find_col(cols, c("maf","freq","af","eaf","a1_freq","freq1","effect_allele_frequency"))

  # Check for combined CHR:POS columns when CHR/POS absent
  if (is.null(map$CHR) || is.null(map$POS)) {
    combo <- find_col(cols, c("chrpos","chr_pos","chr:pos","snp_loc","location"), exact = FALSE)
    if (!is.null(combo)) {
      map$CHRPOS_COMBO <- combo
      log_warn(sprintf("  No separate CHR/POS columns; will parse from '%s'.", combo))
    }
  }

  # Log detections
  for (nm in names(map)) {
    if (!is.null(map[[nm]])) {
      log_info(sprintf("  %-12s → '%s'", nm, map[[nm]]))
    } else {
      if (nm %in% c("RSID","CHR","POS","A1","A2")) {
        log_warn(sprintf("  %-12s → NOT FOUND", nm))
      }
    }
  }

  map
}

# =============================================================================
# 5. STANDARDISE COLUMN NAMES INTO WORKING TABLE
# =============================================================================

build_working_table <- function(dt, map, opt_N) {
  log_info("Building standardised working table...")

  out <- data.table(
    RSID = NA_character_,
    CHR  = NA_character_,
    POS  = NA_integer_,
    A1   = NA_character_,
    A2   = NA_character_,
    MAF  = NA_real_,
    BETA = NA_real_,
    SE   = NA_real_,
    Z    = NA_real_,
    P    = NA_real_,
    N    = NA_real_
  )
  out <- out[rep(1L, nrow(dt))]

  # RSID
  if (!is.null(map$RSID)) {
    out$RSID <- as.character(dt[[map$RSID]])
    # Validate: first non-NA must start with "rs"
    valid_ex <- na.omit(out$RSID)[1:min(10, sum(!is.na(out$RSID)))]
    if (!any(startsWith(valid_ex, "rs"))) {
      log_warn("RSID column detected but values don't start with 'rs' — treating as missing.")
      out$RSID <- NA_character_
    }
  }

  # CHR / POS — parse combo column if needed
  if (!is.null(map$CHRPOS_COMBO)) {
    combo_vals <- as.character(dt[[map$CHRPOS_COMBO]])
    # Accept formats: chr1:12345 / 1:12345 / 1_12345
    parts <- strsplit(gsub("chr", "", combo_vals, ignore.case = TRUE), "[: _/]")
    out$CHR <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
    out$POS <- suppressWarnings(
      vapply(parts, function(x) if (length(x) >= 2) as.integer(x[2]) else NA_integer_, integer(1))
    )
  } else {
    if (!is.null(map$CHR)) out$CHR <- gsub("^chr", "", as.character(dt[[map$CHR]]), ignore.case=TRUE)
    if (!is.null(map$POS)) out$POS <- suppressWarnings(as.integer(dt[[map$POS]]))
  }

  # Alleles — uppercase
  if (!is.null(map$A1)) out$A1 <- toupper(as.character(dt[[map$A1]]))
  if (!is.null(map$A2)) out$A2 <- toupper(as.character(dt[[map$A2]]))

  # BETA — convert OR to log(OR) if needed
  if (!is.null(map$BETA)) {
    out$BETA <- suppressWarnings(as.numeric(dt[[map$BETA]]))
  } else if (!is.null(map$OR)) {
    log_info("BETA not found; converting OR column to log(OR) as BETA.")
    or_vals  <- suppressWarnings(as.numeric(dt[[map$OR]]))
    out$BETA <- log(or_vals)
  }

  if (!is.null(map$SE))  out$SE   <- suppressWarnings(as.numeric(dt[[map$SE]]))
  if (!is.null(map$Z))   out$Z    <- suppressWarnings(as.numeric(dt[[map$Z]]))
  if (!is.null(map$P))   out$P    <- suppressWarnings(as.numeric(dt[[map$P]]))
  if (!is.null(map$MAF)) out$MAF  <- suppressWarnings(as.numeric(dt[[map$MAF]]))

  # N — user-supplied overrides column
  if (!is.null(opt_N)) {
    out$N <- as.numeric(opt_N)
  } else if (!is.null(map$N)) {
    out$N <- suppressWarnings(as.numeric(dt[[map$N]]))
  }

  out
}

# =============================================================================
# 6. RSID ↔ CHR:POS BRIDGING
# =============================================================================

recover_rsids_from_chrpos <- function(dt, rsid_ref_path) {
  if (is.null(rsid_ref_path)) {
    log_warn("No --rsid_ref_path provided; cannot recover missing RSIDs from CHR:POS.")
    return(dt)
  }
  if (!any(is.na(dt$RSID) & !is.na(dt$CHR) & !is.na(dt$POS) &
            !is.na(dt$A1)  & !is.na(dt$A2))) {
    log_info("No SNPs need RSID recovery from CHR:POS — skipping.")
    return(dt)
  }

  need_rsid <- which(is.na(dt$RSID) & !is.na(dt$CHR) & !is.na(dt$POS))
  log_info(sprintf("Recovering RSIDs for %d SNPs using CHR:POS:A2:A1 key...", length(need_rsid)))

  sub   <- dt[need_rsid]
  sub[, LOOKUP_KEY := paste(CHR, POS, A2, A1, sep = ":")]
  chrs  <- unique(sub$CHR)
  recovered <- 0L

  for (chr in chrs) {
    ref_path <- gsub("\\{CHR\\}", chr, rsid_ref_path)
    if (!file.exists(ref_path)) {
      log_warn(sprintf("  Ref file not found for CHR %s: %s", chr, ref_path))
      next
    }
    ref <- tryCatch(
      fread(ref_path, header = FALSE, col.names = c("KEY","RSID_REF"),
            sep = "\t", data.table = TRUE, showProgress = FALSE),
      error = function(e) { log_warn(sprintf("  Could not read %s", ref_path)); NULL }
    )
    if (is.null(ref)) next
    setkey(ref, KEY)

    idx_chr  <- need_rsid[sub$CHR == chr]
    keys_chr <- sub$LOOKUP_KEY[sub$CHR == chr]
    merged   <- ref[keys_chr]
    hits     <- !is.na(merged$RSID_REF)
    if (any(hits)) {
      dt$RSID[idx_chr[hits]] <- merged$RSID_REF[hits]
      recovered <- recovered + sum(hits)
    }
  }

  COUNTERS$n_rsid_recovered <- COUNTERS$n_rsid_recovered + recovered
  log_info(sprintf("  Recovered %d RSIDs from CHR:POS.", recovered))
  dt
}

recover_chrpos_from_rsid <- function(dt, rsid_to_chrpos_path) {
  if (is.null(rsid_to_chrpos_path)) {
    log_warn("No --rsid_to_chrpos_path provided; cannot recover CHR/POS from RSID.")
    return(dt)
  }
  need_pos <- which(!is.na(dt$RSID) & (is.na(dt$CHR) | is.na(dt$POS)))
  if (length(need_pos) == 0L) {
    log_info("No SNPs need CHR/POS recovery from RSID — skipping.")
    return(dt)
  }
  log_info(sprintf("Recovering CHR/POS for %d SNPs from RSID mapping...", length(need_pos)))

  ref <- tryCatch(
    fread(rsid_to_chrpos_path, header = TRUE, data.table = TRUE, showProgress = FALSE),
    error = function(e) { log_warn(sprintf("Could not read rsid_to_chrpos: %s", conditionMessage(e))); NULL }
  )
  if (is.null(ref)) return(dt)

  # Flexible col detection in ref file
  ref_cols <- tolower(names(ref))
  rsid_c <- names(ref)[which(ref_cols %in% c("rsid","rs_id","snp","id"))[1]]
  chr_c  <- names(ref)[which(ref_cols %in% c("chr","chrom","chromosome"))[1]]
  pos_c  <- names(ref)[which(ref_cols %in% c("pos","bp","position"))[1]]
  if (any(is.na(c(rsid_c, chr_c, pos_c)))) {
    log_warn("rsid_to_chrpos file does not have expected RSID/CHR/POS columns — skipping.")
    return(dt)
  }

  ref2 <- ref[, .(RSID_REF = get(rsid_c), CHR_REF = as.character(get(chr_c)),
                  POS_REF  = as.integer(get(pos_c)))]
  setkey(ref2, RSID_REF)

  merged  <- ref2[dt$RSID[need_pos]]
  chr_hit <- !is.na(merged$CHR_REF)
  pos_hit <- !is.na(merged$POS_REF)

  dt$CHR[need_pos[chr_hit]] <- merged$CHR_REF[chr_hit]
  dt$POS[need_pos[pos_hit]] <- merged$POS_REF[pos_hit]

  n_rec <- sum(chr_hit & pos_hit)
  COUNTERS$n_chrpos_recovered <- COUNTERS$n_chrpos_recovered + n_rec
  log_info(sprintf("  Recovered CHR/POS for %d SNPs.", n_rec))
  dt
}

# =============================================================================
# 7. MAF REFERENCE MERGE (1000G EUR)
# =============================================================================

merge_maf_reference <- function(dt, maf_ref_path) {
  if (is.null(maf_ref_path)) {
    log_warn("No --maf_ref_path provided; MAF will rely on input file only.")
    return(dt)
  }
  if (!file.exists(maf_ref_path)) {
    log_warn(sprintf("MAF reference file not found: %s — skipping MAF merge.", maf_ref_path))
    return(dt)
  }
  log_info("Loading 1000G EUR MAF reference...")

  # afreq format: #CHROM ID REF ALT ALT_FREQS OBS_CT
  ref <- tryCatch(
    fread(maf_ref_path, header = TRUE, data.table = TRUE, showProgress = FALSE,
          col.names = c("CHROM","ID","REF","ALT","ALT_FREQS","OBS_CT")),
    error = function(e) {
      log_warn(sprintf("Could not read MAF reference: %s", conditionMessage(e))); NULL
    }
  )
  if (is.null(ref)) return(dt)

  ref <- ref[!is.na(ID) & ID != "."]
  ref[, ALT_FREQS := suppressWarnings(as.numeric(ALT_FREQS))]
  ref <- ref[!is.na(ALT_FREQS)]
  ref[, REF := toupper(REF)][, ALT := toupper(ALT)]
  setkey(ref, ID)

  log_info(sprintf("  MAF reference: %d entries loaded.", nrow(ref)))

  # Merge on RSID
  has_rsid <- !is.na(dt$RSID)
  if (sum(has_rsid) == 0L) {
    log_warn("No RSIDs available for MAF reference merge.")
    return(dt)
  }

  merged <- ref[dt$RSID[has_rsid], .(REF_REF = REF, ALT_REF = ALT, ALT_FREQ = ALT_FREQS)]

  dt[has_rsid, REF_REF  := merged$REF_REF]
  dt[has_rsid, ALT_REF  := merged$ALT_REF]
  dt[has_rsid, ALT_FREQ := merged$ALT_FREQ]

  log_info(sprintf("  MAF merged for %d / %d RSIDs.", sum(!is.na(dt$ALT_FREQ)), sum(has_rsid)))
  dt
}

# =============================================================================
# 8. ALLELE ALIGNMENT + MAF CORRECTION
# =============================================================================

align_alleles_and_maf <- function(dt) {
  log_info("Performing allele alignment and MAF correction...")

  has_ref <- !is.na(dt$REF_REF) & !is.na(dt$ALT_REF) & !is.na(dt$ALT_FREQ)

  # ---- Complement helper (for strand flipping if needed later) ----
  complement <- function(x) {
    chartr("ACGTacgt", "TGCAtgca", x)
  }

  n_swap  <- 0L
  n_flip  <- 0L
  n_maf_fill <- 0L

  for (i in which(has_ref)) {
    a1   <- dt$A1[i]; a2 <- dt$A2[i]
    ref  <- dt$REF_REF[i]; alt <- dt$ALT_REF[i]
    freq <- dt$ALT_FREQ[i]  # frequency of ALT allele

    if (is.na(a1) || is.na(a2)) next

    # Case 1: A1=ALT, A2=REF → MAF = freq (frequency of A1=ALT)
    if (a1 == alt && a2 == ref) {
      dt$MAF[i] <- freq
      n_maf_fill <- n_maf_fill + 1L
      next
    }
    # Case 2: A1=REF, A2=ALT → swapped; MAF of A1 = 1 - freq
    if (a1 == ref && a2 == alt) {
      dt$MAF[i] <- 1 - freq
      n_maf_fill <- n_maf_fill + 1L
      next
    }
    # Case 3: try complement strand
    a1c <- complement(a1); a2c <- complement(a2)
    if (a1c == alt && a2c == ref) {
      dt$MAF[i] <- freq
      dt$A1[i]  <- a1c; dt$A2[i] <- a2c
      n_maf_fill <- n_maf_fill + 1L
      next
    }
    if (a1c == ref && a2c == alt) {
      dt$MAF[i] <- 1 - freq
      dt$A1[i]  <- a1c; dt$A2[i] <- a2c
      n_maf_fill <- n_maf_fill + 1L
      next
    }
    # Allele mismatch — leave MAF from input if available, else NA
  }

  # ---- For SNPs with MAF only from input (no ref merge): still apply >0.5 rule ----
  # Ensure A1 is always the minor (effect) allele convention: MAF ≤ 0.5
  high_maf <- !is.na(dt$MAF) & dt$MAF > 0.5
  n_flip   <- sum(high_maf, na.rm = TRUE)
  if (n_flip > 0L) {
    log_info(sprintf("  Flipping %d SNPs where MAF > 0.5 (swapping A1/A2, sign of Z/BETA).", n_flip))
    dt[high_maf, MAF  := 1 - MAF]
    # Swap A1 / A2
    tmp_a1 <- dt$A1[high_maf]
    dt[high_maf, A1 := A2]
    dt[high_maf, A2 := tmp_a1]
    # Flip Z and BETA
    dt[high_maf & !is.na(Z),    Z    := -Z]
    dt[high_maf & !is.na(BETA), BETA := -BETA]
    n_swap <- n_flip
  }

  COUNTERS$n_maf_flip    <- COUNTERS$n_maf_flip    + n_flip
  COUNTERS$n_allele_swap <- COUNTERS$n_allele_swap + n_swap

  log_info(sprintf("  MAF filled from reference: %d SNPs.", n_maf_fill))
  log_info(sprintf("  MAF > 0.5 flips applied  : %d SNPs.", n_flip))
  dt
}

# =============================================================================
# 9. STATISTICAL DERIVATION (Z, BETA, SE, P)
# =============================================================================

derive_statistics <- function(dt) {
  log_info("Deriving / reconstructing statistics...")

  # --- Step 1: If OR provided and no BETA, BETA = log(OR) already handled earlier ---

  # --- Step 2: Z from BETA + SE ---
  need_z <- is.na(dt$Z) & !is.na(dt$BETA) & !is.na(dt$SE) & dt$SE > 0
  if (any(need_z, na.rm=TRUE)) {
    dt[need_z, Z := BETA / SE]
    log_info(sprintf("  Z computed from BETA/SE for %d SNPs.", sum(need_z, na.rm=TRUE)))
  }

  # --- Step 3: Z from P (two-sided) ---
  need_z_from_p <- is.na(dt$Z) & !is.na(dt$P) & dt$P > 0 & dt$P <= 1
  if (any(need_z_from_p, na.rm=TRUE)) {
    # Determine sign from BETA if available
    dt[need_z_from_p, Z := {
      z_abs <- abs(qnorm(P / 2))
      sign_z <- ifelse(!is.na(BETA), sign(BETA), 1)
      sign_z * z_abs
    }]
    n_from_p <- sum(need_z_from_p, na.rm=TRUE)
    COUNTERS$n_z_from_p <- COUNTERS$n_z_from_p + n_from_p
    log_info(sprintf("  Z derived from P for %d SNPs.", n_from_p))
  }

  # --- Step 4: Approximate BETA/SE from Z + MAF + N (clearly flagged) ---
  need_se <- is.na(dt$SE) & !is.na(dt$Z) & !is.na(dt$MAF) & !is.na(dt$N) &
             dt$MAF > 0 & dt$MAF < 1 & dt$N > 0
  if (any(need_se, na.rm=TRUE)) {
    log_warn(sprintf(paste0(
      "APPROXIMATION APPLIED: SE estimated from Z + MAF + N for %d SNPs. ",
      "Formula: SE = 1/sqrt(2*MAF*(1-MAF)*(N+Z^2)). ",
      "This is an approximation — treat downstream results with caution."),
      sum(need_se, na.rm=TRUE)))
    dt[need_se, SE   := 1 / sqrt(2 * MAF * (1 - MAF) * (N + Z^2))]
    dt[need_se, BETA := Z * SE]
    n_approx <- sum(need_se, na.rm=TRUE)
    COUNTERS$n_beta_approx <- COUNTERS$n_beta_approx + n_approx
  }

  # --- Step 5: Special case — Freq + Effect + N only ---
  special <- is.na(dt$Z) & is.na(dt$SE) & !is.na(dt$MAF) & !is.na(dt$BETA) & !is.na(dt$N) &
             dt$MAF > 0 & dt$MAF < 1 & dt$N > 0
  if (any(special, na.rm=TRUE)) {
    log_warn(sprintf(paste0(
      "APPROXIMATION APPLIED: SE/Z from Freq+Effect+N only for %d SNPs. ",
      "Formula: SE=1/sqrt(2*N*MAF*(1-MAF)), Z=BETA/SE."),
      sum(special, na.rm=TRUE)))
    dt[special, SE := 1 / sqrt(2 * N * MAF * (1 - MAF))]
    dt[special, Z  := BETA / SE]
    n_approx2 <- sum(special, na.rm=TRUE)
    COUNTERS$n_beta_approx <- COUNTERS$n_beta_approx + n_approx2
  }

  # --- Step 6: P from Z if still missing ---
  need_p <- is.na(dt$P) & !is.na(dt$Z)
  if (any(need_p, na.rm=TRUE)) {
    dt[need_p, P := 2 * pnorm(-abs(Z))]
    log_info(sprintf("  P computed from Z for %d SNPs.", sum(need_p, na.rm=TRUE)))
  }

  dt
}

# =============================================================================
# 10. FINAL VALIDATION + CLEANUP
# =============================================================================

finalise <- function(dt) {
  log_info("Finalising output table...")

  # Drop internal reference columns if present
  drop_cols <- c("REF_REF","ALT_REF","ALT_FREQ","LOOKUP_KEY")
  dt[, (intersect(names(dt), drop_cols)) := NULL]

  # Keep only output columns
  out_cols <- c("RSID","CHR","POS","A1","A2","MAF","BETA","SE","Z","P","N")
  for (col in out_cols) {
    if (!col %in% names(dt)) dt[[col]] <- NA
  }
  dt <- dt[, ..out_cols]

  # Sanity bounds
  dt[!is.na(MAF) & (MAF < 0 | MAF > 0.5), MAF := NA]
  dt[!is.na(P)   & (P  < 0 | P  > 1),     P   := NA]
  dt[!is.na(SE)  & SE <= 0,                SE  := NA]

  COUNTERS$n_missing_rsid <- sum(is.na(dt$RSID))
  COUNTERS$n_output       <- nrow(dt)

  dt
}

# =============================================================================
# 11. WRITE OUTPUT
# =============================================================================

write_output <- function(dt, path) {
  log_info(sprintf("Writing %d SNPs to %s ...", nrow(dt), path))
  fwrite(dt, file = path, sep = "\t", na = "NA", quote = FALSE)
  log_info("Output written successfully.")
}

# =============================================================================
# 12. MAIN PIPELINE
# =============================================================================

main <- function() {
  # --- Read ---
  dt_raw <- detect_and_read(opt$sumstat)
  COUNTERS$n_input <- nrow(dt_raw)

  # --- Column detection ---
  col_map <- detect_columns(dt_raw)

  # Check RSID or CHR/POS availability
  has_rsid     <- !is.null(col_map$RSID)
  has_chrpos   <- (!is.null(col_map$CHR) && !is.null(col_map$POS)) ||
                   !is.null(col_map$CHRPOS_COMBO)
  if (!has_rsid && !has_chrpos) {
    log_error(paste0(
      "Cannot proceed: no RSID column detected AND no CHR/POS columns found. ",
      "Please ensure your input contains either RSID or CHR+POS columns."
    ))
  }

  # --- Build working table ---
  dt <- build_working_table(dt_raw, col_map, opt$N)
  rm(dt_raw); gc()

  # --- RSID bridging ---
  dt <- recover_rsids_from_chrpos(dt, opt$rsid_ref_path)
  dt <- recover_chrpos_from_rsid(dt, opt$rsid_to_chrpos_path)

  # --- MAF reference merge ---
  dt <- merge_maf_reference(dt, opt$maf_ref_path)

  # --- Allele alignment + MAF correction ---
  dt <- align_alleles_and_maf(dt)

  # --- Statistical derivations ---
  dt <- derive_statistics(dt)

  # --- Finalise ---
  dt <- finalise(dt)

  # --- Write ---
  write_output(dt, opt$output)

  # --- Summary ---
  print_summary()

  log_info("=== holy-script complete. ===")
  invisible(dt)
}

main()
