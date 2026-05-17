#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# easy-GWAS/lib/run_gwas.sh — 单类型双软件 GWAS
# ═══════════════════════════════════════════════════════════════
set -euo pipefail

VCF="$1"; PHENO_CSV="$2"; TYPE="$3"; OUTDIR="$4"; TRAIT_COL="$5"
GCTA="$6"; GEMMA="$7"
MAF_SNP="$8"; GENO_SNP="$9"; MAF_INDEL="${10}"; GENO_INDEL="${11}"; MAF_SV="${12}"; GENO_SV="${13}"

case "$TYPE" in
    SNP|snp)   MAF=$MAF_SNP; GENO=$GENO_SNP; LABEL="SNP" ;;
    INDEL|indel) MAF=$MAF_INDEL; GENO=$GENO_INDEL; LABEL="INDEL" ;;
    SV|sv)     MAF=$MAF_SV; GENO=$GENO_SV; LABEL="SV" ;;
    *) echo "ERROR: type must be SNP/INDEL/SV"; exit 1 ;;
esac

TDIR="$OUTDIR/$LABEL"
mkdir -p "$TDIR"/{filtered,plink,gcta,gemma}
LOG="$TDIR/${LABEL}.log"
exec > >(tee -a "$LOG") 2>&1

echo "══════════════════════════════════════════════"
echo "  easy-GWAS — $LABEL  |  $(date)"
echo "  MAF>$MAF  missing<$GENO  trait_col=$TRAIT_COL"
echo "══════════════════════════════════════════════"

ERRORS=0

# ══ 1. Filter VCF ══
echo "[1/5] Filtering VCF..."
BCFTOOLS="${BCFTOOLS:-$(which bcftools 2>/dev/null || echo bcftools)}"
PLINK="${PLINK:-$(which plink2 2>/dev/null || echo plink2)}"
FILT_VCF="$TDIR/filtered/${LABEL}.filtered.vcf.gz"
$BCFTOOLS view -i "F_MISSING < $GENO" -q ${MAF}:minor -m2 -M2 \
    -Oz -o "$FILT_VCF" "$VCF" --threads 4
$BCFTOOLS index -f "$FILT_VCF"

# ══ 2. VCF → PLINK ══
echo "[2/5] Converting to PLINK..."
if [ "$TYPE" = "SNP" ]; then
    ID_FMT="@:#:\$r:\$a"
else
    ID_FMT="@:#"
fi
$PLINK --vcf "$FILT_VCF" --set-all-var-ids $ID_FMT \
    --allow-extra-chr \
    --make-bed --out "$TDIR/plink/${LABEL}_tmp" --threads 4 --silent

# Sort samples to match phenotype order
python3 - "$PHENO_CSV" "$TDIR" << 'PYEOF'
import sys, csv
pheno_csv, tdir = sys.argv[1], sys.argv[2]
with open(pheno_csv) as f:
    reader = csv.reader(f)
    next(reader)
    ids = [row[0] for row in reader]
with open(f"{tdir}/sample_order.txt", 'w') as f:
    for iid in ids:
        f.write(f"0 {iid}\n")
PYEOF

# Check overlap before sorting
N_PHENO=$(wc -l < "$TDIR/sample_order.txt")
N_VCF=$(wc -l < "$TDIR/plink/${LABEL}_tmp.fam")
VCF_IDS=$(awk '{print $2}' "$TDIR/plink/${LABEL}_tmp.fam")
PHENO_IDS=$(awk '{print $2}' "$TDIR/sample_order.txt")
MATCHED=$(comm -12 <(echo "$VCF_IDS" | sort) <(echo "$PHENO_IDS" | sort) | wc -l)
echo "  VCF: $N_VCF samples, Pheno: $N_PHENO samples, Overlap: $MATCHED"
if [ "$MATCHED" -eq 0 ]; then
    echo "ERROR: No samples overlap between VCF and phenotype!"
    echo "VCF first 5 IDs: $(echo "$VCF_IDS" | head -5 | tr '\n' ' ')"
    echo "Pheno first 5 IDs: $(echo "$PHENO_IDS" | head -5 | tr '\n' ' ')"
    exit 1
fi

$PLINK --bfile "$TDIR/plink/${LABEL}_tmp" \
    --keep "$TDIR/sample_order.txt" \
    --indiv-sort f "$TDIR/sample_order.txt" \
    --allow-extra-chr \
    --make-bed --out "$TDIR/plink/${LABEL}" --threads 4 --silent

rm -f "$TDIR"/plink/${LABEL}_tmp.* "$TDIR"/sample_order.txt
N_VAR=$(wc -l < "$TDIR/plink/${LABEL}.bim")
N_SAM=$(wc -l < "$TDIR/plink/${LABEL}.fam")
echo "  $N_VAR variants × $N_SAM samples"

# ══ 3. Phenotype ══
echo "[3/5] Preparing phenotype..."
python3 - "$PHENO_CSV" "$TDIR" "$LABEL" "$TRAIT_COL" << 'PYEOF'
import sys, csv
pheno_csv, tdir, label, trait_col = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

fam_ids = []
with open(f"{tdir}/plink/{label}.fam") as f:
    for l in f:
        fam_ids.append(l.strip().split()[1])

with open(pheno_csv) as f:
    reader = csv.reader(f)
    next(reader)
    phe_dict = {}
    for row in reader:
        if len(row) > trait_col and row[trait_col] != 'NA':
            phe_dict[row[0]] = row[trait_col]

# GCTA format: FID IID trait (header)
with open(f"{tdir}/gcta/pheno.txt", 'w') as f:
    f.write(f"FID IID {label}\n")
    for iid in fam_ids:
        f.write(f"0 {iid} {phe_dict.get(iid, '-9')}\n")

# GEMMA format: single column, matched to FAM order
with open(f"{tdir}/gemma/pheno_single.txt", 'w') as f:
    for iid in fam_ids:
        f.write(f"{phe_dict.get(iid, 'NA')}\n")

n_valid = sum(1 for v in phe_dict.values())
print(f"  {n_valid} phenotyped, {len(fam_ids)} total")
PYEOF

# ══ 4. GCTA LOCO GWAS ══
echo "[4/5] GCTA LOCO GWAS..."
if $GCTA --mlma-loco --bfile "$TDIR/plink/${LABEL}" \
    --pheno "$TDIR/gcta/pheno.txt" \
    --out "$TDIR/gcta/gwas" --threads 8 2>"$TDIR/gcta/gcta.err"; then
    N_HITS=$(awk 'NR>1 && $9<1e-5' "$TDIR/gcta/gwas.loco.mlma" 2>/dev/null | wc -l)
    echo "  GCTA: done ($N_HITS hits at p<1e-5)"
else
    ERRORS=$((ERRORS+1))
    echo "  GCTA: FAILED (see $TDIR/gcta/gcta.err)"
fi

# ══ 5. GEMMA LMM GWAS ══
echo "[5/5] GEMMA GWAS..."
awk '{$6=0; print}' OFS='\t' "$TDIR/plink/${LABEL}.fam" > "$TDIR/gemma/${LABEL}.fam"
cp "$TDIR/plink/${LABEL}.bed" "$TDIR/gemma/${LABEL}.bed"
cp "$TDIR/plink/${LABEL}.bim" "$TDIR/gemma/${LABEL}.bim"

cd "$TDIR/gemma"
if $GEMMA -bfile "${LABEL}" -gk 1 -o "${LABEL}_kin" 2>gemma_kin.err; then
    GEMMA_OUT=$($GEMMA -bfile "${LABEL}" -k "output/${LABEL}_kin.cXX.txt" \
        -lmm 1 -miss 1.0 -maf 0 -p pheno_single.txt -o "${LABEL}_gwas" 2>gemma_gwas.err) || true
    PVE=$(echo "$GEMMA_OUT" | grep "pve estimate" | awk -F'=' '{gsub(/ /,""); print $2}')
    N_HITS=$(awk 'NR>1 && $12<1e-5' "output/${LABEL}_gwas.assoc.txt" 2>/dev/null | wc -l)
    echo "  GEMMA: done, pve=$PVE ($N_HITS hits at p<1e-5)"
    # Move output up
    cp "output/${LABEL}_gwas.assoc.txt" "${LABEL}_gwas.assoc.txt" 2>/dev/null || true
else
    ERRORS=$((ERRORS+1))
    echo "  GEMMA: FAILED (see $TDIR/gemma/gemma_*.err)"
fi
cd - > /dev/null

# ══ Summary ══
echo ""
echo "═══ $LABEL GWAS complete ═══"
echo "  GCTA:  $TDIR/gcta/gwas.loco.mlma"
echo "  GEMMA: $TDIR/gemma/${LABEL}_gwas.assoc.txt"
if [ $ERRORS -gt 0 ]; then
    echo "  WARNING: $ERRORS tool(s) failed — check *.err files"
fi
echo "Done: $(date)"
