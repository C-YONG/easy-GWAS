#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# easy-GWAS/lib/compute_gec.R — GEC via eigenvalue decomposition (Li & Ji)
#
# Computes effective number of independent tests per chromosome
# using PLINK LD matrix + poolr::meff()
# Output: single integer = total Neff across all chromosomes
# ═══════════════════════════════════════════════════════════════
suppressMessages(library(poolr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript compute_gec.R <plink_prefix> <n_chr> <outfile>")
}

PLINK   <- args[1]   # plink2 binary
BFILE   <- args[2]   # PLINK bed prefix
N_CHR   <- as.integer(args[3])
OUTFILE <- args[4]

total_neff <- 0
chr_neff <- c()

for (chr in 1:N_CHR) {
  ld_file <- paste0("_gec_chr", chr, ".ld.gz")
  
  # PLINK: compute pairwise LD within 5Mb window, r²>0
  cmd <- sprintf("%s --bfile %s --chr %d --r2-unphased --ld-window-kb 5000 --ld-window-r2 0 --out _gec_chr%d --silent 2>/dev/null",
                 PLINK, BFILE, chr, chr)
  system(cmd, ignore.stderr = TRUE)
  
  ld_path <- paste0("_gec_chr", chr, ".vcor")
  if (!file.exists(ld_path)) next
  
  # Read LD pairs
  ld <- tryCatch(read.table(ld_path, header = TRUE), error = function(e) NULL)
  if (is.null(ld) || nrow(ld) < 2) {
    # No LD pairs → all independent
    # Count SNPs on this chromosome from bim
    bim <- read.table(paste0(BFILE, ".bim"), header = FALSE)
    n_snps <- sum(bim$V1 == chr)
    chr_neff[chr] <- n_snps
    total_neff <- total_neff + n_snps
    next
  }
  
  # Get unique SNP list and build correlation matrix
  snps <- unique(c(ld$ID_A, ld$ID_B))
  n <- length(snps)
  
  if (n > 5000) {
    # Too many SNPs: use Li & Ji approximation = sum(sqrt(max(1-λ², 0)))
    # Quick approximation via LD pruning count
    bim <- read.table(paste0(BFILE, ".bim"), header = FALSE)
    n_snps <- sum(bim$V1 == chr)
    # Conservative: use LD pruned count from earlier step
    prune_file <- paste0(dirname(BFILE), "/", basename(BFILE), "_ld.prune.in")
    if (file.exists(prune_file)) {
      prune_lines <- readLines(prune_file)
      prune_ids <- sapply(strsplit(prune_lines, ":"), function(x) x[1])
      n_eff <- sum(as.numeric(prune_ids) == chr)
    } else {
      n_eff <- n_snps
    }
    chr_neff[chr] <- n_eff
    total_neff <- total_neff + n_eff
    next
  }
  
  # Build correlation matrix (sparse → full for small n)
  r_mat <- diag(1, n)
  rownames(r_mat) <- snps
  colnames(r_mat) <- snps
  
  for (i in 1:nrow(ld)) {
    a <- ld$ID_A[i]
    b <- ld$ID_B[i]
    r2 <- ld$UNPHASED_R2[i]
    if (a %in% snps && b %in% snps) {
      r_mat[a, b] <- sqrt(max(r2, 0))
      r_mat[b, a] <- r_mat[a, b]
    }
  }
  
  # Li & Ji method: Neff from eigenvalues
  neff <- tryCatch({
    meff(r_mat, method = "li4j")
  }, error = function(e) {
    # Fallback: count eigenvalues > 1
    eig <- eigen(r_mat, only.values = TRUE)$values
    sum(eig)
  })
  
  chr_neff[chr] <- neff
  total_neff <- total_neff + neff
  
  # Cleanup
  unlink(ld_path)
}

cat(total_neff, "\n")
writeLines(as.character(round(total_neff)), OUTFILE)
cat("Per-chromosome:", paste(round(chr_neff), collapse = " "), "\n")
