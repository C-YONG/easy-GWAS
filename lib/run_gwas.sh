#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# easy-GWAS/lib/run_gwas.sh — 单类型四软件 GWAS
# ═══════════════════════════════════════════════════════════════
set -euo pipefail

VCF="$1"; PHENO_CSV="$2"; TYPE="$3"; OUTDIR="$4"; TRAIT_COL="$5"
GCTA="$6"; GEMMA="$7"; LDAK="$8"; N_PCA="${9:-5}"; METHOD="${10:-all}"; DO_GEC="${11:-0}"
MAF_SNP="${12}"; GENO_SNP="${13}"; MAF_INDEL="${14}"; GENO_INDEL="${15}"; MAF_SV="${16}"; GENO_SV="${17}"

# Helper: check if a method should run
run_method() { [[ "$METHOD" == "all" || ",$METHOD," == *",$1,"* ]]; }

case "$TYPE" in
    SNP|snp)   MAF=$MAF_SNP; GENO=$GENO_SNP; LABEL="SNP" ;;
    INDEL|indel) MAF=$MAF_INDEL; GENO=$GENO_INDEL; LABEL="INDEL" ;;
    SV|sv)     MAF=$MAF_SV; GENO=$GENO_SV; LABEL="SV" ;;
    *) echo "ERROR: type must be SNP/INDEL/SV"; exit 1 ;;
esac

TDIR="$OUTDIR/$LABEL"
mkdir -p "$TDIR"/{filtered,plink,gcta,gemma,ldak}
LOG="$TDIR/${LABEL}.log"
exec > >(tee -a "$LOG") 2>&1

echo "══════════════════════════════════════════════"
echo "  easy-GWAS — $LABEL  |  $(date)"
echo "  MAF>$MAF  missing<$GENO  trait_col=$TRAIT_COL  PCA=$N_PCA"
echo "  methods: $METHOD  gec: $([ "$DO_GEC" = "1" ] && echo on || echo off)"
echo "══════════════════════════════════════════════"

ERRORS=0

# ══ 1. Filter VCF ══
echo "[1/7] Filtering VCF..."
BCFTOOLS="${BCFTOOLS:-$(which bcftools 2>/dev/null || echo bcftools)}"
PLINK="${PLINK:-$(which plink2 2>/dev/null || echo plink2)}"
FILT_VCF="$TDIR/filtered/${LABEL}.filtered.vcf.gz"
$BCFTOOLS view -i "F_MISSING < $GENO" -q ${MAF}:minor -m2 -M2 \
    -Oz -o "$FILT_VCF" "$VCF" --threads 4
$BCFTOOLS index -f "$FILT_VCF"

# ══ 2. VCF → PLINK ══
echo "[2/7] Converting to PLINK..."
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

# Check overlap
N_PHENO=$(wc -l < "$TDIR/sample_order.txt")
N_VCF=$(wc -l < "$TDIR/plink/${LABEL}_tmp.fam")
VCF_IDS=$(awk '{print $2}' "$TDIR/plink/${LABEL}_tmp.fam")
PHENO_IDS=$(awk '{print $2}' "$TDIR/sample_order.txt")
MATCHED=$(comm -12 <(echo "$VCF_IDS" | sort) <(echo "$PHENO_IDS" | sort) | wc -l)
echo "  VCF: $N_VCF samples, Pheno: $N_PHENO samples, Overlap: $MATCHED"
if [ "$MATCHED" -eq 0 ]; then
    echo "ERROR: No samples overlap between VCF and phenotype!"
    exit 1
fi

$PLINK --bfile "$TDIR/plink/${LABEL}_tmp" \
    --keep "$TDIR/sample_order.txt" --indiv-sort f "$TDIR/sample_order.txt" \
    --allow-extra-chr --make-bed --out "$TDIR/plink/${LABEL}" --threads 4 --silent

rm -f "$TDIR"/plink/${LABEL}_tmp.* "$TDIR"/sample_order.txt
N_VAR=$(wc -l < "$TDIR/plink/${LABEL}.bim")
N_SAM=$(wc -l < "$TDIR/plink/${LABEL}.fam")
echo "  $N_VAR variants × $N_SAM samples"

# ══ 3. GEC ══
echo "[3/7] Computing GEC..."
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
if [ "$DO_GEC" = "1" ]; then
    JAVA="$(find "$SCRIPT_DIR/../bin/jre" -name java -type f 2>/dev/null | head -1)"
    if [ -n "$JAVA" ] && [ -f "$SCRIPT_DIR/../bin/gec/gec.jar" ]; then
        $JAVA -jar -Xmx2g "$SCRIPT_DIR/../bin/gec/gec.jar" --effect-number \
            --plink-binary "$TDIR/plink/${LABEL}" --genome --maf 0 \
            --out "$TDIR/plink/gec_out" 2>/dev/null || true
        if [ -f "$TDIR/plink/gec_out.sum" ]; then
            N_EFF=$(awk 'NR==2{print $2}' "$TDIR/plink/gec_out.sum")
            [ -z "$N_EFF" ] && N_EFF="$N_VAR"
        else N_EFF="$N_VAR"; fi
    else
        $PLINK --bfile "$TDIR/plink/${LABEL}" --indep-pairwise 50 5 0.2 \
            --out "$TDIR/plink/${LABEL}_ld" --silent 2>/dev/null || true
        N_EFF=$(wc -l < "$TDIR/plink/${LABEL}_ld.prune.in" 2>/dev/null || echo "$N_VAR")
    fi
    echo "  Neff = $N_EFF (GEC)"
else
    N_EFF="$N_VAR"
    echo "  Neff = $N_EFF (all markers, no GEC)"
fi
echo "$N_EFF" > "$TDIR/plink/.neff"

# ══ 4. Phenotype ══
echo "[4/7] Preparing phenotype..."
python3 - "$PHENO_CSV" "$TDIR" "$LABEL" "$TRAIT_COL" << 'PYEOF'
import sys, csv
pheno_csv, tdir, label, trait_col = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
fam_ids = []
with open(f"{tdir}/plink/{label}.fam") as f:
    for l in f: fam_ids.append(l.strip().split()[1])
with open(pheno_csv) as f:
    reader = csv.reader(f); next(reader)
    phe_dict = {}
    for row in reader:
        if len(row) > trait_col and row[trait_col] != 'NA':
            phe_dict[row[0]] = row[trait_col]
with open(f"{tdir}/gcta/pheno.txt", 'w') as f:
    f.write(f"FID IID {label}\n")
    for iid in fam_ids: f.write(f"0 {iid} {phe_dict.get(iid, '-9')}\n")
with open(f"{tdir}/gemma/pheno_single.txt", 'w') as f:
    for iid in fam_ids: f.write(f"{phe_dict.get(iid, 'NA')}\n")
n_valid = sum(1 for v in phe_dict.values())
print(f"  {n_valid} phenotyped, {len(fam_ids)} total")
PYEOF

# ══ 5. GCTA-LOCO ══
if run_method loco; then
echo "[5/7] GCTA LOCO GWAS..."
if $GCTA --mlma-loco --bfile "$TDIR/plink/${LABEL}" \
    --pheno "$TDIR/gcta/pheno.txt" \
    --out "$TDIR/gcta/gwas_loco" --threads 8 2>"$TDIR/gcta/gcta_loco.err"; then
    echo "  GCTA-LOCO: done"
else
    ERRORS=$((ERRORS+1)); echo "  GCTA-LOCO: FAILED"
fi
else
    echo "[5/7] GCTA-LOCO: skipped"
fi

# ══ 6. GCTA-MLMA ══
if run_method mlma; then
echo "[6/7] GCTA MLMA GWAS..."
if $GCTA --bfile "$TDIR/plink/${LABEL}" --make-grm-bin \
    --out "$TDIR/gcta/grm" --threads 8 2>"$TDIR/gcta/gcta_grm.err" && \
   $GCTA --mlma --bfile "$TDIR/plink/${LABEL}" --grm "$TDIR/gcta/grm" \
    --pheno "$TDIR/gcta/pheno.txt" \
    --out "$TDIR/gcta/gwas_mlma" --threads 8 2>"$TDIR/gcta/gcta_mlma.err"; then
    echo "  GCTA-MLMA: done"
else
    ERRORS=$((ERRORS+1)); echo "  GCTA-MLMA: FAILED"
fi
else
    echo "[6/7] GCTA-MLMA: skipped"
fi

# ══ 7. GEMMA ══
if run_method gemma; then
echo "[7/7] GEMMA GWAS..."
awk '{$6=0; print}' OFS='\t' "$TDIR/plink/${LABEL}.fam" > "$TDIR/gemma/${LABEL}.fam"
cp "$TDIR/plink/${LABEL}.bed" "$TDIR/gemma/${LABEL}.bed"
cp "$TDIR/plink/${LABEL}.bim" "$TDIR/gemma/${LABEL}.bim"
cd "$TDIR/gemma"
if $GEMMA -bfile "${LABEL}" -gk 1 -o "${LABEL}_kin" 2>gemma_kin.err; then
    GEMMA_OUT=$($GEMMA -bfile "${LABEL}" -k "output/${LABEL}_kin.cXX.txt" \
        -lmm 1 -miss 1.0 -maf 0 -p pheno_single.txt -o "${LABEL}_gwas" 2>gemma_gwas.err) || true
    PVE=$(echo "$GEMMA_OUT" | grep "pve estimate" | awk -F'=' '{gsub(/ /,""); print $2}')
    echo "  GEMMA: done, pve=$PVE"
    cp "output/${LABEL}_gwas.assoc.txt" "${LABEL}_gwas.assoc.txt" 2>/dev/null || true
else
    ERRORS=$((ERRORS+1)); echo "  GEMMA: FAILED"
fi
cd - > /dev/null
else
    echo "[7/7] GEMMA: skipped"
fi

# ══ 8. LDAK-KVIK + PCA ══
if run_method kvik; then
if [ -n "${LDAK:-}" ] && [ -x "$LDAK" ]; then
    echo "  LDAK-KVIK + ${N_PCA}PCs..."

    # PCA from SNP genotypes (standard for all types)
    PCA_DIR="$OUTDIR/pca"
    mkdir -p "$PCA_DIR"
    if [ ! -f "$PCA_DIR/covar${N_PCA}.txt" ]; then
        # Compute PCA from current PLINK data (shared across all types)
        $PLINK --bfile "$TDIR/plink/${LABEL}" --pca "$N_PCA" \
            --out "$PCA_DIR/pca" --threads 8 --silent 2>/dev/null || true
        awk -v n="$N_PCA" 'NR>1{printf $1" "$2; for(i=3;i<=2+n;i++) printf " "$i; print ""}' \
            "$PCA_DIR/pca.eigenvec" > "$PCA_DIR/covar${N_PCA}.txt"
    fi

    # Prepare LDAK data (truncate INDEL/SV alleles)
    mkdir -p "$TDIR/ldak"
    cp "$TDIR/plink/${LABEL}.bed" "$TDIR/ldak/${LABEL}.bed"
    cp "$TDIR/plink/${LABEL}.fam" "$TDIR/ldak/${LABEL}.fam"
    if [ "$TYPE" != "SNP" ]; then
        awk 'BEGIN{FS=OFS="\t"} {if(NF>=6){a1=substr($5,1,1);a2=substr($6,1,1);
            if(a1==a2){a1="A";a2="T"};$5=a1;$6=a2};print}' \
            "$TDIR/plink/${LABEL}.bim" > "$TDIR/ldak/${LABEL}.bim"
    else
        cp "$TDIR/plink/${LABEL}.bim" "$TDIR/ldak/${LABEL}.bim"
    fi

    cd "$TDIR/ldak"
    if $LDAK --kvik-step1 "${LABEL}_kvik" --bfile "${LABEL}" \
        --pheno "$TDIR/gcta/pheno.txt" --covar "$PCA_DIR/covar${N_PCA}.txt" \
        --allow-many-predictors YES 2>ldak_step1.err; then
        if $LDAK --kvik-step2 "${LABEL}_kvik" --bfile "${LABEL}" \
            --pheno "$TDIR/gcta/pheno.txt" --covar "$PCA_DIR/covar${N_PCA}.txt" 2>ldak_step2.err; then
            echo "  LDAK-KVIK: done"
        else
            ERRORS=$((ERRORS+1)); echo "  LDAK-KVIK step2: FAILED"
        fi
    else
        ERRORS=$((ERRORS+1)); echo "  LDAK-KVIK step1: FAILED"
    fi
    cd - > /dev/null
else
    echo "  LDAK-KVIK: skipped (not installed)"
fi
else
    echo "  LDAK-KVIK: skipped"
fi

# ══ Summary ══
echo ""
echo "═══ $LABEL GWAS complete ═══"
run_method loco && [ -f "$TDIR/gcta/gwas_loco.loco.mlma" ] && \
    echo "  GCTA-LOCO:   $TDIR/gcta/gwas_loco.loco.mlma"
run_method mlma && [ -f "$TDIR/gcta/gwas_mlma.mlma" ] && \
    echo "  GCTA-MLMA:   $TDIR/gcta/gwas_mlma.mlma"
run_method gemma && [ -f "$TDIR/gemma/${LABEL}_gwas.assoc.txt" ] && \
    echo "  GEMMA:       $TDIR/gemma/${LABEL}_gwas.assoc.txt"
run_method kvik && [ -f "$TDIR/ldak/${LABEL}_kvik.step2.assoc" ] && \
    echo "  LDAK-KVIK:   $TDIR/ldak/${LABEL}_kvik.step2.assoc"
if [ $ERRORS -gt 0 ]; then
    echo "  WARNING: $ERRORS tool(s) failed — check *.err files"
fi
echo "Done: $(date)"
