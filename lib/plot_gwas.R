#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# easy-GWAS/lib/plot_gwas.R — Manhattan + QQ plots with GEC support
#
# Usage:
#   Rscript plot_gwas.R <input.mlma|input.assoc.txt> <output_prefix> [--gec] [--format jpg|png]
#
# Input:  GCTA .mlma (Chr SNP bp p) or GEMMA .assoc.txt (chr rs ps p_wald)
# Output: <prefix>_Manhattan.jpg + <prefix>_QQ.jpg
# ═══════════════════════════════════════════════════════════════

suppressMessages(library(CMplot))

# ── Parse args ──
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_gwas.R <input_file> <output_prefix> [--gec] [--format jpg|png]")
}

input_file  <- args[1]
out_prefix  <- args[2]
use_gec     <- "--gec" %in% args
out_format  <- if ("--png" %in% args) "png" else "jpg"
trait_name  <- basename(out_prefix)

# Create output directory and CD into it (CMplot uses setwd)
out_dir <- dirname(out_prefix)
out_name <- basename(out_prefix)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

# ── Read GWAS data ──
if (grepl("\\.assoc\\.txt$", input_file)) {
  # GEMMA format: chr rs ps n_miss allele1 allele0 af beta se logl_H1 l_remle p_wald
  raw <- read.table(input_file, header = TRUE)
  gwas <- data.frame(
    SNP = raw$rs,
    Chromosome = raw$chr,
    Position = raw$ps,
    p = raw$p_wald,
    stringsAsFactors = FALSE
  )
} else {
  # GCTA format: Chr SNP bp A1 A2 Freq b se p
  raw <- read.table(input_file, header = TRUE)
  gwas <- data.frame(
    SNP = raw$SNP,
    Chromosome = raw$Chr,
    Position = raw$bp,
    p = raw$p,
    stringsAsFactors = FALSE
  )
}

# Remove NA p-values
gwas <- gwas[!is.na(gwas$p) & gwas$p > 0, ]
cat("Markers:", nrow(gwas), "\n")

# ── Thresholds ──
N_MARKERS <- nrow(gwas)
bonf_sig  <- 0.05 / N_MARKERS
bonf_sugg <- 1.0 / N_MARKERS

if (use_gec) {
  # GEC: read Neff computed by PLINK --indep-pairwise 50 5 0.2
  # Look for .neff in {gwas_dir}/../../plink/.neff
  neff_file <- file.path(dirname(input_file), "..", "..", "plink", ".neff")
  if (file.exists(neff_file)) {
    N_EFF <- as.numeric(readLines(neff_file)[1])
    cat("GEC Neff (PLINK --indep-pairwise 50 5 0.2):", N_EFF, "\n")
  } else {
    # Try alternative location
    neff_file <- file.path(dirname(input_file), "..", "plink", ".neff")
    if (file.exists(neff_file)) {
      N_EFF <- as.numeric(readLines(neff_file)[1])
      cat("GEC Neff:", N_EFF, "\n")
    } else {
      N_EFF <- N_MARKERS
      cat("GEC Neff file not found, using raw N:", N_EFF, "\n")
    }
  }
  sig_thresh  <- 0.05 / N_EFF
  sugg_thresh <- 1.0 / N_EFF
  threshold_label <- paste0("GEC-sig (", formatC(sig_thresh, format="e", digits=1), ")")
} else {
  sig_thresh  <- bonf_sig
  sugg_thresh <- bonf_sugg
  threshold_label <- paste0("Bonf-sig (", formatC(sig_thresh, format="e", digits=1), ")")
}

cat(sprintf("Significant:  p < %.2e\n", sig_thresh))
cat(sprintf("Suggestive:   p < %.2e\n", sugg_thresh))

# ── Manhattan ──
CMplot(
  gwas,
  plot.type = "m",
  LOG10 = TRUE,
  threshold = c(sig_thresh, sugg_thresh),
  threshold.lty = c(1, 2),
  threshold.lwd = c(1.5, 1),
  threshold.col = c("red", "orange"),
  amplify = TRUE,
  bin.size = 1e6,
  signal.col = c("red", "darkgreen"),
  signal.cex = c(1.2, 1.2),
  signal.pch = c(19, 19),
  file = out_format,
  file.name = paste0(out_name, "_Manhattan"),
  dpi = 300,
  file.output = TRUE,
  verbose = FALSE,
  width = 14,
  height = 8
)

# ── QQ ──
CMplot(
  gwas,
  plot.type = "q",
  box = FALSE,
  file = out_format,
  file.name = paste0(out_name, "_QQ"),
  dpi = 300,
  conf.int = TRUE,
  conf.int.col = NULL,
  threshold.col = "red",
  threshold.lty = 2,
  file.output = TRUE,
  verbose = FALSE
)

cat("Done: ", out_prefix, "_Manhattan.", out_format, " + _QQ.", out_format, "\n", sep="")
