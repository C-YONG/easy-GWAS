# easy-GWAS

**Three-method GWAS (GCTA-LOCO + GCTA-MLMA + GEMMA-LMM) with GEC support.**

Estimates effective number of independent tests using the official GEC software (Li et al. 2012) and generates Manhattan + QQ plots.

## Install

```bash
git clone https://github.com/C-YONG/easy-GWAS.git
cd easy-GWAS
bash install.sh                            # GCTA + GEMMA + GEC + Java JRE
conda install -c bioconda bcftools plink2
```

## Quick start

```bash
# Single variant type (3 methods)
easy-GWAS --out results/ single snp.vcf.gz pheno.csv SNP

# All three types
easy-GWAS --out results/ batch snp.vcf.gz indel.vcf.gz sv.vcf.gz pheno.csv

# Manhattan + QQ plots (with GEC thresholds)
easy-GWAS plot --gec results/SNP/gcta/gwas_loco.loco.mlma SNP_trait
```

## GWAS methods

| Method | Software | Description |
|--------|----------|-------------|
| LOCO | GCTA | Leave-one-chromosome-out (gold standard) |
| MLMA | GCTA | Mixed linear model with full GRM (faster) |
| LMM | GEMMA | Mixed linear model with Wald test + pve |

## GEC

The official GEC v0.2 (Li et al. 2012, pmglab.top) computes the effective number of independent tests via eigenvalue decomposition of the LD matrix. Used for significance thresholds in plots.

## Output

```
results/SV/
├── gcta/
│   ├── gwas_loco.loco.mlma   GCTA LOCO (Chr bp p)
│   └── gwas_mlma.mlma        GCTA MLMA
├── gemma/
│   └── SV_gwas.assoc.txt     GEMMA LMM (chr ps p_wald)
└── plink/
    └── gec_out.sum           GEC: Neff + thresholds
```

## Filtering defaults

| Type  | MAF  | Missingness |
|-------|------|-------------|
| SNP   | 0.05 | 0.05        |
| INDEL | 0.05 | 0.05        |
| SV    | 0.01 | 0.20        |
