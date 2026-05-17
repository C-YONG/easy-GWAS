# easy-GWAS

**Two-tool GWAS (GCTA-LOCO + GEMMA-LMM) for SNP, INDEL, SV.**

## Install

```bash
git clone https://github.com/C-YONG/easy-GWAS.git
cd easy-GWAS
bash install.sh                          # Downloads GCTA & GEMMA
conda install -c bioconda bcftools plink2
```

## Usage

```bash
# Single variant type
easy-GWAS --out results/ single snp.vcf.gz pheno.csv SNP

# All three types
easy-GWAS --out results/ batch snp.vcf.gz indel.vcf.gz sv.vcf.gz pheno.csv

# Options
--trait N    Trait column in CSV (0=FID, 1=first trait, default: 1)
```

## Output

```
results/SNP/
├── gcta/gwas.loco.mlma      GCTA LOCO results (Chr bp p)
└── gemma/SNP_gwas.assoc.txt GEMMA LMM results (chr ps p_wald)
```

## Filtering defaults

| Type  | MAF  | Missingness |
|-------|------|-------------|
| SNP   | 0.05 | 0.05        |
| INDEL | 0.05 | 0.05        |
| SV    | 0.01 | 0.20        |
