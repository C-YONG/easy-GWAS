# easy-GWAS

**Four-method GWAS (GCTA-LOCO + GCTA-MLMA + GEMMA-LMM + LDAK-KVIK) with GEC support.**

## Install

```bash
git clone https://github.com/C-YONG/easy-GWAS.git
cd easy-GWAS
bash install.sh                            # GCTA + GEMMA + LDAK + GEC + Java JRE
conda install -c bioconda bcftools plink2
```

## Quick start

```bash
# Single variant type (4 methods)
easy-GWAS --pca 5 --out results/ single snp.vcf.gz pheno.csv SNP

# Single method only
easy-GWAS --method loco --out results/ single snp.vcf.gz pheno.csv SNP

# Subset of methods
easy-GWAS --method loco,gemma --out results/ single snp.vcf.gz pheno.csv SNP

# All three types
easy-GWAS --pca 5 --out results/ batch snp.vcf.gz indel.vcf.gz sv.vcf.gz pheno.csv

# Manhattan + QQ plots (with GEC thresholds)
easy-GWAS plot --gec results/SNP/gcta/gwas_loco.loco.mlma SNP_TD
```

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
```

## GEC

The official GEC v0.2 (Li et al. 2012, pmglab.top) computes effective independent tests via eigenvalue decomposition of the LD matrix. Used for significance thresholds in plots.

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
└── plink/
    └── gec_out.sum           GEC Neff + thresholds
```

## Filtering defaults

| Type  | MAF  | Missingness |
|-------|------|-------------|
| SNP   | 0.05 | 0.05        |
| INDEL | 0.05 | 0.05        |
| SV    | 0.01 | 0.20        |
