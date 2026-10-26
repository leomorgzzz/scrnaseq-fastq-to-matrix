#!/usr/bin/env bash
# Pseudoalinear los FASTQ y generar la matriz de conteos (barcodes x genes) sin filtrar
# Usar desde la carpeta base del repo, como: bash scripts/02_kb_count.sh

# Si algo falla detiene el script si un comando falla, una variable no esta definida 
# o si falla un paso del pipeline.
set -euo pipefail

# Rutas de entrada (FASTQ y referencia) y de salida.
FASTQ_DIR=data/raw/pbmc_1k_v3_fastqs
REF_DIR=data/reference
OUT_DIR=data/processed/pbmc_1k_v3
mkdir -p "$OUT_DIR"

# kb count: pseudoalinea R2 contra el indice (-i), asigna transcrtos a genes con el t2g (-g), lee barcode
# y UMI de R1 según 10x v3 (-x), usa 8 hilos (-t) y guarda además la matriz en formato h5ad para
# Scanpy (--head). Los FASTQ van por pares (R1 y luego R2) de cada carril, I1 no se usa porque las
# lecturas ya están separadas por muestra. Sin --filter, así la matriz queda sin filtrar.

kb count \
        -i "$REF_DIR/index.idx" \
        -g "$REF_DIR/t2g.txt" \
        -x 10xv3 \
        -o "$OUT_DIR" \
        -t 8 \
        --h5ad \
        "$FASTQ_DIR/pbmc_1k_v3_S1_L001_R1_001.fastq.gz" \
        "$FASTQ_DIR/pbmc_1k_v3_S1_L001_R2_001.fastq.gz" \
        "$FASTQ_DIR/pbmc_1k_v3_S1_L002_R1_001.fastq.gz" \
        "$FASTQ_DIR/pbmc_1k_v3_S1_L002_R2_001.fastq.gz"

