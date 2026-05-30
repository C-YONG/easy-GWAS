# easy-GWAS

**Four-method GWAS (GCTA-LOCO + GCTA-MLMA + GEMMA-LMM + LDAK-KVIK) with auto-plotting and optional GEC.**

## Install

```bash
git clone https://github.com/C-YONG/easy-GWAS.git
cd easy-GWAS
bash install.sh                            # GCTA + GEMMA + LDAK + GEC + Java JRE
conda install -c bioconda bcftools plink2
```

## Quick start

```bash
# Single variant type (4 methods, auto-plot with Bonferroni)
easy-GWAS --pca 5 --out results/ single snp.vcf.gz pheno.csv SNP

# GEC-corrected thresholds (opt-in)
easy-GWAS --gec --pca 5 --out results/ single snp.vcf.gz pheno.csv SNP

# Single method only
easy-GWAS --method loco --out results/ single snp.vcf.gz pheno.csv SNP

# Subset of methods with GEC
easy-GWAS --gec --method loco,gemma --out results/ single snp.vcf.gz pheno.csv SNP

# All three variant types (auto-plot per type)
easy-GWAS --pca 5 --out results/ batch snp.vcf.gz indel.vcf.gz sv.vcf.gz pheno.csv

# Manual plot (when auto-plot fails or custom prefix needed)
easy-GWAS plot --gec results/SNP/gcta/gwas_loco.loco.mlma SNP_TD
```

Plots are auto-generated after GWAS: `{TYPE}_trait{N}_{method}_Manhattan.jpg` + `_QQ.jpg` in each type's result directory.

## GWAS methods

| Method | Software | Description |
|--------|----------|-------------|
| LOCO | GCTA | Leave-one-chromosome-out (gold standard) |
| MLMA | GCTA | Mixed linear model with full GRM |
| LMM | GEMMA | Mixed linear model with Wald test + pve |
| KVIK | LDAK | Elastic net + PCA (default 5 PCs, Nature Genetics 2024) |

## Options

```
--out DIR     Output directory (default: ./easy-gwas-out)
--trait N     Trait column (0=FID, 1=first trait, default: 1)
--pca N       Number of PCs for LDAK-KVIK (default: 5, 0=skip KVIK)
--method M    GWAS methods: all, loco, mlma, gemma, kvik (comma-separated, default: all)
--gec         Use GEC-corrected significance thresholds (off by default)
```

Without `--gec`, Bonferroni correction (0.05 / N_markers) is used for plot thresholds.
With `--gec`, GEC computes effective independent tests for corrected thresholds.

## GEC

The official GEC v0.2 (Li et al. 2012, pmglab.top) computes effective independent tests via eigenvalue decomposition of the LD matrix. Used for significance thresholds in plots when `--gec` is passed.

## Output

```
results/SV/
├── gcta/
│   ├── gwas_loco.loco.mlma   GCTA LOCO (Chr bp p)
│   └── gwas_mlma.mlma        GCTA MLMA
├── gemma/
│   └── SV_gwas.assoc.txt     GEMMA (chr ps p_wald)
├── ldak/
│   └── SV_kvik.step2.assoc   LDAK-KVIK (Wald_P)
├── plink/
│   └── gec_out.sum           GEC Neff + thresholds (only with --gec)
├── SV_trait1_loco_Manhattan.jpg   Auto-generated plot
├── SV_trait1_loco_QQ.jpg          Auto-generated plot
└── ...
```

## Filtering defaults

| Type  | MAF  | Missingness |
|-------|------|-------------|
| SNP   | 0.05 | 0.05        |
| INDEL | 0.05 | 0.05        |
| SV    | 0.01 | 0.20        |
