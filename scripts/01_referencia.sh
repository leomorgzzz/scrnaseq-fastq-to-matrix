#!/usr/bin/env bash
# Descargar el índice prehecho de kallisto (humano) y su tabla t2g.
# Usar desde la carpeta base del repo, como: bash scripts/01_referencia.sh


# Si algo falla detiene el script si un comando falla, una variable no esta definida 
# o si falla un paso del pipeline.
set -euo pipefail

# Localización donde se guarda la referencia, se crea si no existe.
REF_DIR=data/reference
mkdir -p "$REF_DIR"

# Si el índice y el t2g ya existen, no volver a descargarlos.
if [[ -f "$REF_DIR/index.idx" && -f "$REF_DIR/t2g.txt" ]]; then
        echo "El índice ya existe en $REF_DIR; no se descarga de nuevo."
        exit 0
fi

# Descargar el índice humano prehecho (-d), guardando el indice de kmers (-i)
# y la tabla transcrito -> gen (-g).
kb ref -d human \
	-i "$REF_DIR/index.idx" \
	-g "$REF_DIR/t2g.txt"

