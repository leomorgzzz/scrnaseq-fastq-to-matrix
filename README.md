# scRNA-seq: Convertir FASTQ a matriz de conteos

Práctica de procesamiento de datos crudos de single-cell RNA-seq (10x Genomics) con kallisto | bustools.
## 1. Estructura de los datos

### Dataset
- **Fuente:** [10x Genomics, dataset público "1k PBMCs from a Healthy Donor (v3 chemistry)".](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0)
- **Muestra:** Células mononucleares de sangre periférica, aprox. 1,222 células.
- **Química:** Chromium Single Cell 3' v3 (procesador con Cell Ranger 3.0.0).
- **Secuenciador:** Illumina NovaSeq 6000

### Archivos FASTQ

| Archivo | Lectura | Longitud (pb) | Contenido |
|---|---|---|---|
| `..._L001_I1_001.fastq.gz` | I1 | 8 | Índice de muestra, para identificar la biblioteca dentro de una corrida con varias muestras |
| `..._L001_R1_001.fastq.gz` | R1 | 28 | Cell barcode (16 pb) + UMI (12 pb) |
| `..._L001_R2_001.fastq.gz` | R2 | 91 | cDNA, secuencia del transcrito |

**Nomenclatura** (`pbmc_1k_v3_S1_L001_R1_001.fastq.gz`):
- `S1`: Número de muestra de la corrida.
- `L001` / `L002`: Carril de la celda de flujo. La misma biblioteca se secuenció en dos carriles y se combinan al cuantificar.
- `R1` / `R2`: Lectura 1 y Lectura 2 del mismo fragmento. En esta tecnología (10x) no son simétricas, R1 lee barcode y UMI y R2 lee el cDNA.
- `I1`: Lectura del índice de la muestra.
- `001`: Número de fragmento del archivo.

### Formato de una lectura FASTQ

```
@A00228:279:HFWFVDMXX:1:1101:8486:1000 2:N:0:NCATTACT
NACAAAGTCCCCCCCATAATACAGGGGGAGCCACTTGGGCAGGAGGCAGGGAGGGGTCCATTCCCCCTGGTGGGGCTGGTGGGGAGCTGTA
+
#FFFFFFFFFFFFFFF:FFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF
```

1. **Identificador:** Empieza con `@`, contiene la información de la corrida y la posición del cluster.
2. **Secuencia:** Las bases leidas. Según el archivo, es cDNA (R2), barcode y UMI (R1) o índice (I1).
3. **Separador:** `+`
4. **Calidad:** Un carácter por base, codificado como Phred+33.  `Q = código ASCII − 33`, y `Q = −10·log10(P_error)`. Ejemplos: `F` = Q37 (~1 error en 5,000), `:` = Q25, `#` = Q2.

**Desglose del ID de la lectura:**

```
@A00228:279:HFWFVDMXX:1:1101:8486:1000 1:N:0:NCATTACT
```

| Campo | Valor | Significado |
|---|---|---|
| Equipo | `A00228` | Número de serie del secuenciador |
| Corrida | `279` | Número de corrida de ese equipo |
| Flowcell | `HFWFVDMXX` | ID de la celda de flujo |
| Carril | `1` | Carril de la celda de flujo |
| Tile : x : y | `1101:8486:1000` | Sección y coordenadas del cluster en la celda de flujo |
| Lectura | `1` | Número de lectura (`1` en R1, `2` en R2) |
| Filtrada | `N` | `N` = pasó el filtro de calidad; `Y` = fue filtrada |
| Control | `0` | `0` = no es lectura de control |
| Índice | `NCATTACT` | Secuencia del índice de la muestra (la misma que en I1) |

### Estructura del oligo de captura (10x 3' v3)

```
Perla ── [Read 1 primer] ── [Cell barcode 16 pb] ── [UMI 12 pb] ── [poly(dT) ~30 nt]
                             igual en toda la perla   distinto en cada oligo   captura la cola poli-A del mRNA

Biblioteca final:
P5 ─ Read 1 ─ Barcode ─ UMI ─ poly(dT) ─ cDNA ─ Read 2 primer ─ índice i7 ─ P7
      └─ R1: 28 pb ──┘                    └─ R2: 91 pb          └─ I1: 8 pb

```

**Ejemplo con una lectura real de R1:**

```
NGTGATTAGCTGTACT | CGTATGTAAGGT
 barcode (16 pb) |  UMI (12 pb)
```

- **Cell barcode (16 pb):** Identifica la célula. Cada perla tiene un barcode único y cae en una gota con, idealmente una sola célula, así que todas las moléculas capturadas en esa gota comparten barcode.
- **UMI (12 pb):** Identifica la molécula de mRNA original. Permine colapsar los duplicados por PCR, lecturas con el mismo barcode, UMI y gen se cuenta como una sola molécula.
- **poly(dT):** Se une a la cola poli-A de los mRNA. Como la captura empieza en el extremo 3', las lecturas se concentran cerca del 3' del gen.

### Conteo de lecturas

| Carril | Lecturas R1 | Lecturas R2 | Lecturas I1 |
|---|---|---|---|
| L001 | 33,436,697 | 33,436,697 | 33,436,697 |
| L002 | 33,165,190 | 33,165,190 | 33,165,190 |
| **Total** | 66,601,887 | 66,601,887 | 66,601,887 |

- Lecturas por célula, aproximadamente: 66,601,887  / 1,222 células = 54,502 (web summary: 54,502)


### Comandos utilizados

```bash
# Ver el contenido del .tar sin extraerlo
tar -tvf pbmc_1k_v3_fastqs.tar

# Extraer los FASTQ
tar -xvf pbmc_1k_v3_fastqs.tar

# Ver las primeras 2 lecturas (8 líneas) de cada archivo sin descomprimir
zcat pbmc_1k_v3_S1_L001_I1_001.fastq.gz | head -8
zcat pbmc_1k_v3_S1_L001_R1_001.fastq.gz | head -8
zcat pbmc_1k_v3_S1_L001_R2_001.fastq.gz | head -8

# Distribución de longitudes en las primeras 100 lecturas (línea 2 de cada bloque de 4)
zcat pbmc_1k_v3_S1_L001_R1_001.fastq.gz | head -400 | awk 'NR % 4 == 2 {print length($0)}' | sort | uniq -c

# Contar líneas; lecturas = líneas / 4
zcat pbmc_1k_v3_S1_L001_R1_001.fastq.gz | wc -l
```

## 2. Referencia

### Índice utilizado
- **Fuente:** índice prehecho de kallisto descargado con `kb ref -d human` ([pachterlab/kallisto-transcriptome-indices](https://github.com/pachterlab/kallisto-transcriptome-indices)).
- **Versión:** Ensembl 108, ensamble GRCh38; construido con kb-python.
- **Tipo:** *standard*: solo transcritos maduros (cDNA, sin intrones). Incluye genes codificantes, lncRNA y genes de inmunoglobulinas y receptores de células T (IG/TR). La alternativa *nac* agrega transcritos nacientes (con intrones), útil para single-nucleus o RNA velocity.
- **Tamaño:** índice 393 MB, t2g 17 MB.

### Tabla transcrito → gen (t2g)

```
ENST00000343575.11  ENSG00000107562.17  CXCL12  CXCL12-201  10  44377067  44385097  -
ENST00000374426.6   ENSG00000107562.17  CXCL12  CXCL12-202  10  44375867  44385043  -
ENST00000374429.6   ENSG00000107562.17  CXCL12  CXCL12-203  10  44370165  44385092  -
```

| Columna | Contenido |
|---|---|
| 1 | ID de Ensembl del transcrito (`ENST…`; el sufijo `.11` es la versión) |
| 2 | ID de Ensembl del gen (`ENSG…`) |
| 3 | Nombre del gen |
| 4 | Nombre del transcrito (isoforma) |
| 5 | Cromosoma |
| 6 | Posición de inicio |
| 7 | Posición de fin |
| 8 | Hebra (`+` / `-`) |

- **Transcritos:** 227,665
- **Genes:** 39,546
- **Transcritos por gen (promedio):** Aproximadamente 5.8
- **Ejemplo:** CXCL12 tiene 8 isoformas (CXCL12-201 a CXCL12-208), todas del mismo gen `ENSG00000107562`.


### Comandos utilizados

```bash
# Descargar el índice prehecho y el t2g (ver scripts/01_referencia.sh)
bash scripts/01_referencia.sh

# Tamaño de los archivos
ls -lh data/reference/

# Primeras líneas del t2g
head data/reference/t2g.txt

# Número de transcritos (una línea p
wc -l data/reference/t2g.txt

# Número de genes distintos (valores únicos de la columna 2)
cut -f2 data/reference/t2g.txt | sor

# Isoformas de un gen específico
grep -w CXCL12 data/reference/t2g.txt
```
