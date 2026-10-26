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

## 3. Cuantificación con kb count

### Comando

```bash
bash scripts/02_kb_count.sh
```

El script recibe los FASTQ de los dos carriles, el índice de kallisto y el t2g, y genera una matriz de conteos **barcodes × genes sin filtrar**. Los FASTQ se pasan por pares de cada carril (R1 y R2) porque kb espera el barcode y el UMI en el primer archivo de cada par y el cDNA en el segundo. I1 no se usa: el índice de muestra solo sirve para separar bibliotecas de una corrida, y ya vienen separados.

### Pasos internos de kb count

| Paso | Herramienta | Qué hace |
|---|---|---|
| 1 | `kallisto bus` | Pseudoalinea cada R2 contra el índice y guarda, por lectura, su barcode, su UMI y su clase de equivalencia en un archivo BUS (*Barcode, UMI, Set*). |
| 2 | `bustools correct` | Compara cada barcode con la lista de barcodes válidos de 10x v3 y corrige los que difieren por una base (errores de secuenciación). |
| 3 | `bustools sort` | Ordena los registros por barcode y UMI, para que las lecturas de la misma molécula queden juntas. |
| 4 | `bustools count` | Colapsa las lecturas con el mismo barcode + UMI en una sola molécula (elimina duplicados de PCR) y usa el t2g para asignarla a un gen. Genera la matriz. |

**Estructura de lectura para 10x v3** (`kb --list`):

| Elemento | Archivo | Posición |
|---|---|---|
| Cell barcode | R1 (archivo 0) | bases 0–16 |
| UMI | R1 (archivo 0) | bases 16–28 |
| cDNA | R2 (archivo 1) | completa |

### Archivos de salida

| Archivo | Contenido |
|---|---|
| `output.bus` / `output.unfiltered.bus` | Registros binarios de barcode + UMI + clase de equivalencia,antes de contar. Son intermedios y s
| `matrix.ec` | Las clases de equivalencia: qué conjunto de transcritos corresponde a cada número de   clase. |
| `transcripts.txt` | Lista de los transcritos del índice, en el orden que usa `matrix.ec`. |          | `run_info.json` | Resumen del pseuesadas y porcentaje pseudoalineado. |
| `inspect.json` | Resumen de barcodes y UMIs: cuántos hay, cuántos son válidos, medianas por barcode. |
| `counts_unfiltered/cells_x_genes.mtx` | La matriz de conteos en formato disperso (Matrix Market): solo guarda los valores distintos de
| `counts_unfiltered/cells_x_genes.barcodes.txt` | Los barcodes, uno por fila de la matriz. |
| `counts_unfiltered/cells_x_genes.g` | Los genes (IDs de Ensembl y nombres), uno por columna de la matriz. |
| `counts_unfiltered/adata.h5ad` | La misma matriz, con barcodes y genes, en formato AnnData para Scanpy. |

### Pseudoalineamiento (`run_info.json`)

| Métrica | kb (kallisto) | Cell Ranger 3.0.0 |
|---|---|---|
| Lecturas procesadas | 66,601,887 | 66,601,887 |
| Pseudoalineadas / mapeadas al transcriptoma | 38,261,433 (57.4%) | 53.7% |
| Compatibles con un solo transcrito | 14.7% | — |
| Mapeadas a intrones | — (índice sin intrones) | 31.1% |
| Mapeadas a regiones intergénicas | — | 4.8% |


### Barcodes y UMIs (`inspect.json`)

| Métrica | Valor |
|---|---|
| Barcodes distintos (antes de corregir) | 504,233 |
| Barcodes en la lista válida | 230,319 (45.7%) |
| Lecturas con barcode en la lista válida | 95.8% (Cell Ranger: 97.4%) |
| Mediana de lecturas por barcode | 3 |
| Mediana de UMIs por barcode | 1 |
| Barcodes en la matriz sin filtrar | 241,214 |
| Células estimadas por Cell Ranger | 1,222 |
| Mediana de UMIs por célula (Cell Ranger) | 6,628 |


### Conclusión

Se procesaron las 66.6 millones de leudoalineamiento (57.4%) y de barcodes válidos (95.8% de las lecturas) son comparables a las de Cell Ranger. La diferencia principal se
explica por las lecturas intrónicas,incluye. La matriz tiene 241,214 barcodes, en su mayoría gotas vacías; el siguiente paso es separar las células reales con un knee plot.

## 4. Filtrado de células: knee plot

Notebook: `notebooks/01_knee_plot.ipynb`

### Problema

La matriz sin filtrar tiene 241,214 barcodes, pero el experimento capturó alrededor de 1,200 células. La mayoría de los barcodes corresponde a gotas vacías que capturaron RNA ambiental (de células rotas en la suspensión). Hay que decidir qué barcodes son células reales.

### Método

1. Se calcularon los UMIs totales por barcode (suma de cada fila de la matriz).
2. Se ordenaron de mayor a menor y se graficaron contra su rango, en escala log-log (*knee plot*).
3. Se eligió un umbral de UMIs en la zona de caída brusca, que separa las células (muchos UMIs) de las gotas vacías (pocos UMIs).

Se usa escala logarítmica porque los valores abarcan varios órdenes de magnitud: de 1 a alrededor de 59,000 UMIs y de 1 a alrededor de 241,000 barcodes.

![Knee plot](results/figures/knee_plot.png)

La gráfica muestra tres zonas:
- **Meseta superior** (rango 1 a cerca de 1,000): células, con miles a decenas de miles de UMIs.
- **Caída brusca** (rango cercano a 1,150 y 1,350): transición entre células y gotas vacías.
- **Cola inferior** (rango > 1,500): gotas vacías, con 1 a 50 UMIs.

### Elección del umbral

| Rango del barcode | UMIs totales |
|---|---|
| 1,000 | 4,502 |
| 1,150 | 1,612 |
| 1,200 | 805 |
| 1,222 | 661 |
| 1,300 | 230 |
| 1,350 | 116 |
| 1,500 | 53 |

| Umbral (UMIs) | Barcodes retenidos |
|---|---|
| ≥ 100 | 1,369 |
| ≥ 200 | 1,315 |
| **≥ 450** | **1,258** |
| ≥ 500 | 1,251 |
| ≥ 1,000 | 1,177 |
| ≥ 2,000 | 1,131 |

**Umbral elegido: 450 UMIs.** Está en la parte más vertical de la caída, donde se separan claramente las células de las gotas vacías. Entre los rangos 1,150 y 1,350 los UMIs bajan de 1,600 a 116, así que cualquier umbral entre 200 y 1,000 retiene entre 1,177 y 1,315 barcodes. El resultado es robusto a la elección exacta del umbral.

### Resultado

| | Valor |
|---|---|
| Barcodes antes del filtro | 241,214 |
| Células retenidas | 1,258 |
| Barcodes eliminados | 239,956 (99.5%) |
| Células estimadas por Cell Ranger 3.0.0 | 1,222 |
| Diferencia | +36 (+2.9%) |

La matriz filtrada se guardó en `data/processed/pbmc_1k_v3/adata_filtered.h5ad`.

