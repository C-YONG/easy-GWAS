#!/bin/bash
# easy-GWAS/install.sh — 下载 GCTA + GEMMA + GEC + JRE
set -euo pipefail
BINDIR="$(cd "$(dirname "$0")" && pwd)/bin"
mkdir -p "$BINDIR"
cd "$BINDIR"

echo "═══ Installing easy-GWAS dependencies ═══"

if [ ! -f gcta64 ]; then
    echo "  Downloading GCTA..."
    curl -sLO "https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-4-x86_64.zip"
    unzip -qo gcta-1.94.1-linux-kernel-4-x86_64.zip
    GCTA_BIN=$(find gcta-1.94.1-linux-kernel-4-x86_64 -type f -name 'gcta*' ! -name '*.zip' | head -1)
    cp "$GCTA_BIN" gcta64 && chmod +x gcta64
    rm -rf gcta-1.94.1* __MACOSX
    echo "  GCTA: $(./gcta64 2>&1 | head -1)"
fi

if [ ! -f gemma ]; then
    echo "  Downloading GEMMA..."
    curl -sLO "https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz"
    gunzip -f gemma-0.98.5-linux-static-AMD64.gz
    mv gemma-0.98.5-linux-static-AMD64 gemma && chmod +x gemma
    echo "  GEMMA: $(./gemma 2>&1 | head -1)"
fi

if [ ! -f gec/gec.jar ]; then
    echo "  Downloading GEC..."
    curl -sLO "https://pmglab.top/gec/data/archive/v0.2/gecV0.2.zip"
    unzip -qo gecV0.2.zip && rm gecV0.2.zip
    if [ ! -d jre ]; then
        echo "  Downloading Java JRE (~41MB)..."
        curl -sLO "https://github.com/adoptium/temurin11-binaries/releases/download/jdk-11.0.25%2B9/OpenJDK11U-jre_x64_linux_hotspot_11.0.25_9.tar.gz"
        mkdir -p jre && cd jre
        tar xzf ../OpenJDK11U-jre_x64_linux_hotspot_11.0.25_9.tar.gz
        rm ../OpenJDK11U-jre_x64_linux_hotspot_11.0.25_9.tar.gz
        cd ..
    fi
    JAVA=$(find jre -name java -type f | head -1)
    echo "  GEC: $($JAVA -jar gec/gec.jar 2>&1 | grep "GEC!" | head -1)"
fi

echo ""; echo "  Note: bcftools + plink2: conda install -c bioconda bcftools plink2"
echo "═══ Done ═══"
