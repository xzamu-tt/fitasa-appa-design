# Reporte técnico — Diseño computacional de variantes termoestables de la fitasa AppA

**Laboratorio de Fisicoquímica e Ingeniería de Proteínas**
**PI:** Daniel Alejandro Fernández
**Autor del trabajo:** Xzamu
**Fecha del reporte:** 2026-04-16
**Período de ejecución:** Etapa I (jun 2024 – ene 2025): análisis estructural y de conservación · Etapa II (jun – jul 2025): rediseño generativo, validación y selección final.
**Alcance:** Documento exclusivamente sobre la fase **computacional** del proyecto.

---

## 1. Resumen

Se diseñaron *in silico* variantes de la **fitasa AppA de _Escherichia coli_** (PDB **1DKL**, clase HAP) con el objetivo de aumentar su termoestabilidad sin comprometer la arquitectura del sitio activo ni la afinidad por fitato. El pipeline combinó: (i) análisis de conservación evolutiva con **ConSurf** sobre 172 homólogos únicos, (ii) delimitación del sitio activo a **7 Å del ligando** (sobre estructura análoga 1DKP), (iii) construcción de un set de **~201 posiciones redesignables** (fijando sitio activo + residuos más conservados, con objetivo de **~70 % de identidad de secuencia**, siguiendo la lógica de Sumida et al., 2024), (iv) rediseño con **ProteinMPNN** a tres temperaturas (0.1, 0.2, 0.3) con sesgo negativo a cisteína, (v) predicción estructural de las 51 secuencias con **AlphaFold 3 Server**, (vi) relajación y scoring físico con **PyRosetta FastRelax / ref2015**, y (vii) ranking multimétrica combinada pLDDT+Rosetta en dos escenarios (50/50 y 60/40). Se seleccionaron **3 variantes finales** para síntesis génica: **Fitasa-01, Fitasa-02 y Fitasa-03**.

---

## 2. Contexto y objetivo

### 2.1 Objetivo

Diseñar variantes de la fitasa AppA de _E. coli_ con **mayor termoestabilidad** sin afectar su **afinidad catalítica por el sustrato** (fitato / InsP6), utilizando herramientas de IA de inverse folding y filtros físico-químicos. Se seleccionan **3 variantes candidatas** para expresión y purificación (fase experimental posterior, fuera del alcance de este reporte).

### 2.2 Relevancia

Las fitasas industriales deben sobrevivir las temperaturas de peletizado en alimento balanceado (**~65–95 °C** por decenas de segundos). AppA combina alta actividad catalítica con termoestabilidad limitada. La redistribución de residuos no esenciales en superficie/core, preservando sitio activo y motivos conservados, es una estrategia validada para mejorar la tolerancia térmica (Sumida et al., 2024; Navone et al., 2021; Xing et al., 2023).

### 2.3 Inspiración metodológica

Todo el pipeline se basa en **Sumida KH et al. (2024). _Improving Protein Expression, Stability, and Function with ProteinMPNN_. JACS 146(3):2054–2061**, ajustado a un target fitasa. Interpretación clave del artículo: el objetivo **no es fijar el 70 % de los residuos**, sino alcanzar un **~70 % de identidad de secuencia final** respecto al WT. Como ProteinMPNN recupera espontáneamente una fracción de los residuos no fijados (*sequence recovery* ≈ 48–52 % en backbones nativos), se puede **fijar menos del 70 %** y aun así quedar cerca de esa identidad.

---

## 3. Pipeline computacional

```
┌────────────────────────┐
│ PDB 1DKL (AppA, E.coli)│
└──────────┬─────────────┘
           │
           ▼
┌────────────────────────────────────┐
│ (1) Conservación (ConSurf)         │  172 homólogos · score Bayesiano
│ (2) Sitio activo (7 Å de ligando)  │  1DKP como referencia estructural
│ (3) Definición de set redesignable │  ~201 residuos (70 % identidad target)
└──────────┬─────────────────────────┘
           │
           ▼
┌────────────────────────────────────┐
│ (4) Rediseño con ProteinMPNN       │  3 temperaturas (0.1 / 0.2 / 0.3)
│     Modelo v_48_020 · bias Cys-3.0 │  17 seqs/temp = 51 totales
└──────────┬─────────────────────────┘
           │
           ▼
┌────────────────────────────────────┐
│ (5) Plegamiento con AlphaFold 3    │  50 estructuras predichas
│     pLDDT promedio por secuencia   │
└──────────┬─────────────────────────┘
           │
           ▼
┌────────────────────────────────────┐
│ (6) Threading + FastRelax Rosetta  │  ref2015 · top-10 por pLDDT
│     Energía total por modelo       │
└──────────┬─────────────────────────┘
           │
           ▼
┌────────────────────────────────────┐
│ (7) Ranking combinado (2 pesos)    │  50/50 y 60/40 pLDDT–Rosetta
│     Selección de 3 finales         │
└──────────┬─────────────────────────┘
           │
           ▼
┌────────────────────────────────────┐
│ Fitasa-01, Fitasa-02, Fitasa-03    │  → síntesis génica
└────────────────────────────────────┘
```

---

## 4. Paso 1 — Estructura base

- **PDB seleccionado:** [1DKL](https://www.rcsb.org/structure/1DKL) (AppA de _E. coli_).
- **Archivo:** `Escritorio/Protein_Engineering/1dkl.pdb` (copia local) y `Fitasa_PyRosetta/input/fitasa_backbone_clean_1DKL.pdb` (limpio para Rosetta).
- **Uso:** tal cual del PDB, sin minimización previa ni reconstrucción de residuos faltantes.
- **Plegamiento:** dos dominios (α/β + α) con bolsillo catalítico positivamente cargado entre ellos.
- **Motivos catalíticos:** **RHGXRXP** (His nucleófilo) y **HD/HAE** (donador protónico).
- **Extensión:** ~410 residuos por cadena; 1DKL es homodímero, la cadena A fue la usada como base.

**Candidatos descartados:**

| PDB   | Motivo de descarte                                |
|-------|---------------------------------------------------|
| 4TSR  | Mutante, referencia principal no publicada al momento |
| 1DKQ  | Mutante con sustitución en la His catalítica → sitio activo alterado |

Se usó **1DKP** únicamente como referencia para definir el sitio activo (cocristal con análogo de fitato), **no** como backbone de diseño.

---

## 5. Paso 2 — Análisis de conservación (ConSurf)

### 5.1 Servidor y parámetros

**Servidor:** [ConSurf-DB](https://consurfdb.tau.ac.il/main_output_new.php?pdb=1DKL&chain=A)
**Cita:** Ashkenazy H et al. _ConSurf 2016_. NAR, 2016 · Ben Chorin A et al. _ConSurf-DB_. Prot. Sci. 2020.

| Parámetro                        | Valor                                  |
|----------------------------------|----------------------------------------|
| Estructura de entrada            | 1DKL, cadena A                         |
| Base de datos                    | UNIREF90                               |
| Búsqueda de homólogos            | HMMER (E-value 1e-4, 1 iteración)      |
| CD-HIT cutoff                    | 95 % (identidad máxima entre homólogos)|
| Coverage mínimo                  | 60 % de la query                       |
| MSA                              | MAFFT                                  |
| Filogenia                        | Neighbor Joining con distancias ML     |
| Modelo de sustitución            | WAG (mejor ajuste)                     |
| Scoring de conservación          | Método Bayesiano                       |

### 5.2 Resultados

- **1 934** homólogos recuperados → **180** pasan filtros → **172 únicos CD-HIT** usados en el cálculo.
- Distancia pareada promedio: **1.27** (rango 0.05 – 2.05).
- Archivo base: `1DKL_B_msa_aa_variety_percentage.csv` — porcentaje de conservación por posición.

### 5.3 Derivados

- Script **`ConsVisual.py`** (carpeta `MSA-CONSENSUS/GitHub_ProteinConsensusSequence/Visualizador-Conservacion/`) para mapear conservación al **B-factor** y colorear en PyMOL (espectro `blue_white_red`).
- Sesiones PyMOL generadas:
  - `CONSURF_ConservedAA/consurf_1DKL_A/1DKL_consurf_pymol_session.pse`
  - `CONSURF_ConservedAA/consurf_1DKL_A/1DKL_consurf_CBS_pymol_session.pse`
  - `MSA-CONSENSUS/Conservados.pse`
- Análisis específico sobre **prolinas, cisteínas y puentes disulfuro** (objetivo: decidir si se permitía su modificación dado su peso estructural).

---

## 6. Paso 3 — Sitio activo (7 Å del ligando)

### 6.1 Criterio

Usando **1DKP** (complejo fitasa–ligando análogo de fitato) como referencia, se extrajeron todos los residuos con al menos un átomo a **≤ 7 Å del ligando**. Estos residuos nunca se redistribuyen.

### 6.2 Residuos identificados (39 posiciones)

Archivo: `FIXED_RESIDUES/Active_Site_1DKP/Active_site_7aFromLigand.txt`

```
ARG16  ALA17  GLY18  ARG20  ALA21  PRO22  THR23  LYS24
GLY45  ASP88  ASP90  ARG92  THR93  PHE125 ASN126 LYS129
ASN204 SER212 LEU213 SER215 MET216 LEU217 GLU219 ILE220
LEU223 HIS250 GLN253 PHE254 GLN258 ARG267 GLY302 HIS303
ASP304 THR305 ASN306 ASP325 THR327 PRO328 PRO329
```

Incluye los motivos canónicos de HAP: el núcleo **R16-H17-G18-R20** (parte de RHGXRXP) y **H303-D304** (HD).

### 6.3 Sesión PyMOL

`FIXED_RESIDUES/Active_Site_1DKP/1DKP_ActiveSite.pse`

---

## 7. Paso 4 — Construcción del set redesignable

### 7.1 Lógica de selección (principio 70 % identidad)

El objetivo no era fijar literalmente el 70 %, sino garantizar que la **secuencia diseñada final compartiera ~70 % de identidad** con el WT. Dado que ProteinMPNN recupera espontáneamente ~48–52 % de los residuos no fijados (overall_confidence y seq_rec reportados por la herramienta al final de cada diseño), el set fijado puede ser menor al 70 % y aun así alcanzar esa identidad.

### 7.2 Procedimiento

1. Ordenar todas las posiciones por conservación (ConSurf).
2. Tomar las **50 más conservadas** del ranking.
3. **Excluir** de ese Top 50 las posiciones que ya forman parte del **sitio activo** (para no contarlas dos veces): 10 posiciones **(A21, A126, A129, A204, A212, A213, A217, A253, A258, A267)**. Archivo: `Final_FixedAA/ActiveSiteAA-ExcludedFromTop50.csv`.
4. **Unir**: sitio activo (39) + Top-50 conservadas no-activo (~40) + otras posiciones conservadas complementarias hasta llegar al umbral deseado, **dejando ~201 posiciones libres para rediseño**.
5. Exportar la lista de posiciones **redesignables** al formato requerido por ProteinMPNN: `FIXED_RESIDUES/Final_FixedAA/FIXED_residues_format.txt` (pese al nombre del archivo, esta lista contiene las posiciones **a rediseñar**, es decir el complemento al set fijo).

### 7.3 Resultado

- **Posiciones rediseñables:** 201 (desde A1 hasta A410).
- **Posiciones fijas:** ~209 (sitio activo + residuos más conservados).
- **Identidad esperada:** ~70 % entre variante diseñada y WT, tras recuperación espontánea de MPNN.

### 7.4 Dos escenarios paralelos: BN vs FLLS

Se corrió el pipeline con **dos sets de residuos** en paralelo, guardados en:

- `ProteinMPNN-outputs_BN/` — set "**BN**" (bien / set final)
- `PoteinMPNN-outputs_FLLS/` — set "**FLLS**" (fallos / set que no cumplió criterios)

> Los argumentos exactos del script (`ArgumentosProteinMPNN`) son **idénticos** en ambos; la diferencia real estuvo en el contenido del archivo `fixed_residues` cargado en ejecución. La rama **BN fue la que produjo los candidatos finales**; la rama FLLS se descartó por peores resultados aguas abajo (pLDDT / Rosetta más bajos).

---

## 8. Paso 5 — Rediseño con ProteinMPNN

### 8.1 Configuración (script `ArgumentosProteinMPNN`)

```bash
python run.py \
    --model_type "protein_mpnn" \
    --seed 111 \
    --pdb_path "./inputs/1dkl.pdb" \
    --temperature {0.1 | 0.2 | 0.3} \
    --out_folder "./outputs/redesign_residuesT0{1|2|3}" \
    --redesigned_residues "A1 A2 A4 A8 ... A410"  # 201 posiciones
    --bias_AA "C:-3.0" \
    --batch_size 2 \
    --number_of_batches 8
```

- **Modelo:** ProteinMPNN checkpoint **`proteinmpnn_v_48_020.pt`** (0.20 Å de ruido de entrenamiento, version 48).
  *Nota:* el header de salida reporta `use_ligand_context=True`; este flag es informativo del wrapper común con LigandMPNN, pero el modelo usado es **ProteinMPNN** (no LigandMPNN). No se incorporó explícitamente el fitato en el diseño.
- **Seed fijo (111)** → reproducibilidad.
- **Bias `C:-3.0`**: penaliza fuertemente introducir cisteínas nuevas (evita formación de disulfuros artefactuales que puedan alterar plegamiento).
- **Sampling:** 3 temperaturas (T01 = 0.1 conservador, T02 = 0.2 balanceado, T03 = 0.3 diverso).
- **Producción:** `batch_size=2 × number_of_batches=8 = 16` secuencias por temperatura + referencia WT → **17 secuencias/temperatura**.
- **Total:** **51 secuencias diseñadas** (17 × 3 temperaturas).

### 8.2 Métricas reportadas por MPNN

Cada secuencia se acompaña de:
- `overall_confidence` ≈ 0.40–0.50 (probabilidad conjunta media bajo el modelo)
- `ligand_confidence` ≈ igual a overall_confidence (sin ligando activo)
- `seq_rec` ≈ **0.45–0.52** (fracción de residuos recuperados respecto al WT en las posiciones redesignables), consistente con lo reportado por Dauparas et al. 2022.

### 8.3 Outputs

- `ProteinMPNN-outputs_BN/redesign_residuesT0X/seqs/1dkl.fa` (FASTA por temperatura).
- `ProteinMPNN-outputs_BN/redesign_residuesT0X/backbones/` (backbones re-empaquetados por MPNN).

---

## 9. Paso 6 — Plegamiento con AlphaFold 3

### 9.1 Configuración

| Parámetro     | Valor                                              |
|---------------|----------------------------------------------------|
| Servidor      | [AlphaFold Server](https://alphafoldserver.com/)   |
| Input         | Secuencia de proteína únicamente (sin ligando, sin MSA de usuario) |
| Modelos       | Configuración por defecto (5 semillas, output JSON + CIF) |
| Output        | `AlphaFold-T01-T03/` — 50 subcarpetas (1 por secuencia) |

### 9.2 Cálculo de pLDDT promedio

Script R: `ProteinMPNN-outputs_BN/AlphaFold3-Analysis.R`

- Se leen los archivos `full_data_<n>.json` de cada subcarpeta.
- Se extrae `atom_plddts` y se calcula el **promedio sobre todos los átomos** del modelo.
- Resultados ordenados descendentemente.

### 9.3 Resultados agregados

Archivo: `resultados_promedios_plddt.csv` (50 modelos).

**Rango general de pLDDT:** 94.93 – 96.19 (todos **≫ 85**, umbral de Sumida et al.; todos los diseños son de alta confianza estructural).

**Distribución por temperatura (pLDDT promedio):**

| Temperatura | Secuencias | pLDDT medio | Comentario                       |
|-------------|-----------:|------------:|----------------------------------|
| T01 (0.1)   |         17 |      ~95.7  | Los más conservadores            |
| T02 (0.2)   |         17 |      ~95.6  | Balanceados                      |
| T03 (0.3)   |         17 |      ~95.8  | Diversos pero plegables          |

**Top 10 por pLDDT** (archivo `top_10_secuencias_plddt.csv`):

| Rank | Modelo                       | Grupo | pLDDT   |
|-----:|------------------------------|-------|--------:|
|  1   | t03_2025_07_03_11_59_6       | T03   | 96.1949 |
|  2   | t01_2025_07_03_11_45         | T01   | 96.1948 |
|  3   | t03_2025_07_03_11_59_3 ★     | T03   | 96.1779 |
|  4   | t02_2025_07_03_11_48 ★       | T02   | 96.0333 |
|  5   | t02_2025_07_03_11_51_2       | T02   | 96.0194 |
|  6   | t03_2025_07_03_12_00_8       | T03   | 95.9607 |
|  7   | t03_2025_07_03_11_59_4       | T03   | 95.8787 |
|  8   | t02_2025_07_03_11_48_2       | T02   | 95.8746 |
|  9   | t02_2025_07_03_11_50_3 ★     | T02   | 95.8485 |
|  10  | t02_2025_07_03_11_49_3       | T02   | 95.8480 |

★ = seleccionadas finalmente para síntesis.

---

## 10. Paso 7 — Relajación y scoring físico con PyRosetta

### 10.1 Protocolo (`scritps/evaluate_sequence.py`)

1. **Backbone plantilla:** `fitasa_backbone_clean_1DKL.pdb` (WT limpiado).
2. **Threading:** se monta cada una de las Top-10 secuencias sobre el backbone WT con `MutateResidue`.
3. **Relajación:** `FastRelax()` con `scorefxn = ref2015` e init flags `-ex1 -ex2 -use_input_sc` (rotámeros extendidos).
4. **Output:** PDB individual por modelo en `output/threaded_structures/` + score total en `output/scores.csv`.
5. **Múltiples modelos por diseño:** cada secuencia AF3 produce 5 semillas ⇒ Rosetta las evalúa todas y se usa **ranking_top1.tsv** para quedarse con la mejor pose por diseño.

### 10.2 Cita

- Alford RF et al. _The Rosetta all-atom energy function (ref2015)_. JCTC 2017.
- Chaudhury S et al. _PyRosetta_. PLoS ONE 2010.

### 10.3 Resultados top Rosetta (mejor modelo por secuencia)

| Rank | Diseño                       | Modelo                            | Rosetta E (REU) |
|-----:|------------------------------|------------------------------------|----------------:|
|  1   | t02_2025_07_03_11_50_3 ★     | t02_2025_07_03_11_50_3_3           |       -812.332  |
|  2   | t02_2025_07_03_11_48 ★       | t02_2025_07_03_11_48_0             |       -803.833  |
|  3   | t02_2025_07_03_11_51_2       | t02_2025_07_03_11_51_2_4           |       -795.720  |
|  4   | t02_2025_07_03_11_50_3 ★     | t02_2025_07_03_11_50_3_2           |       -791.613  |
|  7   | t03_2025_07_03_11_59_3 ★     | t03_2025_07_03_11_59_3_0           |       -784.871  |
|  ..  |                              |                                    |                 |

Se aprecia que los **tres candidatos finales** ocupan lugares consistentes en el ranking Rosetta (energías entre −784 y −812 REU, claramente más favorables que el WT relajado como referencia interna).

---

## 11. Paso 8 — Ranking combinado pLDDT + Rosetta

Script R embebido en `AlphaFold3-Analysis.R` (sección "ANÁLISIS DE DATOS FINAL").

### 11.1 Normalización

- `plddt_scaled`: (pLDDT − min) / (max − min) — 1 = mejor.
- `rosetta_scaled`: (max_E − E) / (max_E − min_E) — 1 = más negativo = mejor.

### 11.2 Dos escenarios de ponderación

| Escenario | w_pLDDT | w_Rosetta | Preferencia                         |
|-----------|--------:|----------:|-------------------------------------|
| 50/50     |    0.50 |      0.50 | Balanceado                          |
| 60/40     |    0.40 |      0.60 | Prioriza energía Rosetta            |

### 11.3 Selección final: **intersección/consenso**

Los 3 candidatos sintetizados son los que aparecen de forma **consistente** en los primeros rangos de **ambos** escenarios (Top-12 de 50/50 y Top-12 de 60/40), privilegiando adicionalmente la diversidad de temperaturas (no los tres con la misma T).

**Archivos asociados:**
- `Fitasa_PyRosetta/input/top12_por_escenario.csv`
- `Fitasa_PyRosetta/input/top50_por_escenario.csv`
- `Fitasa_PyRosetta/input/rankings_top6_por_escenario.csv`

---

## 12. Variantes finales seleccionadas

Archivo de secuencias: `ProteinMPNN-outputs_BN/SecuenciasFinales_Fitasa.fasta`.

| ID lab     | Timestamp MPNN                 | T de diseño | pLDDT (AF3) | Rosetta E (REU, mejor pose) | Notas                                              |
|------------|--------------------------------|-------------|------------:|---------------------------:|----------------------------------------------------|
| **Fitasa-01** | `t03_2025_07_03_11_59_3` | 0.3         |     96.1779 |                    −784.871 | Máxima diversidad (T03) y buen consenso pLDDT+Rosetta |
| **Fitasa-02** | `t02_2025_07_03_11_48`   | 0.2         |     96.0333 |                    −803.833 | Alta pLDDT + 2ª mejor energía Rosetta              |
| **Fitasa-03** | `t02_2025_07_03_11_50_3` | 0.2         |     95.8485 |                    −812.332 | **Mejor energía Rosetta** del set                  |

### 12.1 Secuencias completas

Ver Anexo A (al final del documento) y `SecuenciasFinales_Fitasa.fasta`.

### 12.2 Perfil de mutaciones vs WT

Dado el objetivo de ~70 % de identidad, cada variante contiene aproximadamente **120–130 sustituciones** respecto al WT AppA (1DKL). El detalle posición a posición puede extraerse con:
```bash
python -c "from Bio import pairwise2; ..."  # alineamiento pairwise
```
o desde el CSV de identidades que se genere a partir de `SecuenciasFinales_Fitasa.fasta`.

---

## 13. Análisis paralelos de soporte (no entran al filtro final pero justifican el pipeline)

| Carpeta                              | Contenido                                                  | Uso                                                                                       |
|--------------------------------------|------------------------------------------------------------|-------------------------------------------------------------------------------------------|
| `HAPYs-Comparison/`                  | MSA, sitios activos y LigPlot de 12 fitasas HAP            | Justificar conservación del motivo RHGXRXP y definición del sitio de unión                |
| `PTPhys-Comparison/`                 | Comparación con fosfatasas tipo PTP                        | Contextualizar la diferencia mecanística entre clases de fitasas                           |
| `ESTRUCTURA-HAPS/`                   | Artículos estructurales + PDB Novel Consensus Bacterial 6-Phytase | Antecedente de diseño de consenso                                                       |
| `MSA-CONSENSUS/GitHub_ProteinConsensusSequence/` | Script `ConsVisual.py` y visualización de conservación  | Mapeo ConSurf → B-factor PyMOL                                                             |
| `CONSURF/`                           | Outputs completos de ConSurf                               | Listas de conservación por residuo                                                         |

---

## 14. Archivos PyMOL (.pse) de referencia

| Archivo                                                        | Contenido                                          |
|----------------------------------------------------------------|----------------------------------------------------|
| `FIXED_RESIDUES/CONSURF_ConservedAA/.../1DKL_consurf_pymol_session.pse`       | Mapa de conservación ConSurf sobre 1DKL            |
| `FIXED_RESIDUES/Active_Site_1DKP/1DKP_ActiveSite.pse`          | Residuos del sitio activo (7 Å del ligando)        |
| `MSA-CONSENSUS/Conservados.pse`                                | Conservación mapeada sobre 1DKL con esferas en P/C |
| `HAPYs-Comparison/SitiosActivos/SESIÓN_COMPARACIÓN_ACTIVOS.pse`| Superposición de sitios activos de 12 HAP          |
| `HAPYs-Comparison/AlrededorDelLigando/AlrededorDelLigando.pse` | Residuos del entorno del ligando en 12 HAP         |
| `HAPYs-Comparison/LigPlot-Comparison/LigPlot-Comparison.pse`   | Interacciones proteína–ligando tipo LigPlot         |

Figuras PNG generadas en `PROYECTO FITASA/pymol_figures/` (cada una con su script `.pml` reproducible; regenerar con `pymol -cq figX.pml`):

| # | Archivo                              | Contenido                                                          |
|---|--------------------------------------|--------------------------------------------------------------------|
| 1 | `fig1_consurf_conservation.png`      | Mapa de conservación ConSurf sobre 1DKL (rojo = conservado, azul = variable) |
| 2 | `fig2_active_site.png`               | Sitio activo a 7 Å del ligando; motivos R16-H17-G18-R20 y H303-D304 resaltados |
| 3 | `fig3_fixed_vs_redesigned.png`       | Cα en esferas: rojo = fijados (~209), azul = rediseñados (201)     |
| 4 | `fig4_wt_vs_variants.png`            | Superposición WT (gris) vs Fitasa-01 (salmón), Fitasa-02 (azul), Fitasa-03 (verde) |

---

## 15. Lo que NO se hizo y quedó como línea futura

A diferencia del plan inicial (archivo `Fitasa.md`), la ronda final **no incluyó**:

| Herramienta             | Propósito                                           | Prioridad recomendada |
|-------------------------|-----------------------------------------------------|-----------------------|
| **FoldX**               | ΔΔG por mutación                                    | Alta                  |
| **Rosetta cartesian_ddg** | ΔΔG físico cartesiano por mutación                 | Alta                  |
| **DynaMut2**            | ΔΔG + modos normales                                | Media                 |
| **GROMACS MD**          | Simulaciones 37 °C y 60 °C, 100 ns, CHARMM36        | Media                 |
| **AutoDock Vina / GNINA** | Docking de fitato sobre variantes                  | Media                 |
| **PLACER**              | Refinamiento del sitio activo por AI                | Baja                  |
| **AF2 ensemble**        | Análisis de rigidez por múltiples semillas          | Baja                  |

Se justifica su omisión por priorizar un filtro **doble** (geometría AF3 + energía Rosetta) suficiente para elegir 3 candidatos dentro del presupuesto de síntesis.

---

## 16. Referencias

1. **Sumida KH et al. 2024.** Improving Protein Expression, Stability, and Function with ProteinMPNN. _JACS_ 146(3):2054–2061. [PMC10811672](https://pmc.ncbi.nlm.nih.gov/articles/PMC10811672/) — **Referencia metodológica primaria**.
2. **Dauparas J et al. 2022.** ProteinMPNN: robust deep learning–based protein sequence design. _Science_. [PMC9997061](https://pmc.ncbi.nlm.nih.gov/articles/PMC9997061/)
3. **Ashkenazy H et al. 2016.** ConSurf 2016. _NAR_.
4. **Ben Chorin A et al. 2020.** ConSurf-DB. _Prot. Sci._ 29:258–267.
5. **Abramson J et al. 2024.** AlphaFold 3. _Nature_.
6. **Alford RF et al. 2017.** The Rosetta ref2015 energy function. _JCTC_. [PMC5717763](https://pmc.ncbi.nlm.nih.gov/articles/PMC5717763/)
7. **Chaudhury S et al. 2010.** PyRosetta. _PLoS ONE_.
8. **Scott BM et al. 2024.** Structural and functional profile of phytases across the domains of life. _Biochimie_. [PMC10982552](https://pmc.ncbi.nlm.nih.gov/articles/PMC10982552/)
9. **Navone L et al. 2021.** Disulfide bond engineering of AppA phytase for increased thermostability. [PMC8010977](https://pmc.ncbi.nlm.nih.gov/articles/PMC8010977/)
10. **Xing H et al. 2023.** Thermostability enhancement of _E. coli_ phytase by directed evolution. [PMC10101328](https://pmc.ncbi.nlm.nih.gov/articles/PMC10101328/)
11. **De Jong JA et al. 2017.** Stability of four commercial phytase products under increasing thermal conditioning temperatures. [PMC7205339](https://pmc.ncbi.nlm.nih.gov/articles/PMC7205339/)

---

## Anexo A — Secuencias de las variantes finales

```
>Fitasa-01 (t03_2025_07_03_11_59_3) | T=0.3 | pLDDT=96.18 | Rosetta=-784.871
SAPALKLEAVVIVSRHGVRAPTKPTPLQRAVTPLPWPEWPVKPGTLTPRGAELIASLGRYWR
ERLVALGLLAASGCPAPGEVQILADVDERTRETGRAFAAGLAPDCDITVHTKPDTTAPDPL
FNPLKTGVCTLDVAAVRDAILAAAGGSIEAFEAARRAEFEALQAVLEFDKSPLCLNRSCDL
TELLPSELVVTPSNVSLTGAVSLASMLSEIFLLQQAQGLPNPGWGRITTPEQWDTLLALHN
AQFYLLQRTPEVARPRATPLLVLIRTALTPHPPRVERYGITLPTKVLYISGHDTNLANLAG
ALELDWSLPGQPDLTPPGGELVFERWRRLADNAYYIQVSAVFQTLEQMRAQTPLSLATPPA
SVPLTLAGCTARNALGYCSLEDFSALVDAAIVPACAL

>Fitasa-02 (t02_2025_07_03_11_48) | T=0.2 | pLDDT=96.03 | Rosetta=-803.833
SEPELVLETVVIVSRHGVRAPTKPTPLMLAVTPHPWPEWPVKPGTLTPRGAQLVAQLGAYWR
ERLVALGLLAAEGCPAPGEVQILADVDERTRRTGEAFAAGLAPDCALPVATRPDTSTPDPL
FNPLKTGVCTLDTAAVTDAILAAAGGSIAAFEAARAAEFDLLEAVLDFEGSPLCLDSACRL
TALLPSTLIVTPTNVSLTGAVSLASMLAEIFLLQQAQGLPDPGWGRITDPAQWDALLSLHN
AQFYLLQRTPEVARPRATPLLELIAAALTPAPPAPQRYGLVLPTRVLYISGHDTNLANLAG
ALELDWELPGQPDRTPPGGELVFERWRRVADNAYWIQVSAVFQTLEQMRNQTPLSLAAPPA
EVPLTLAGCTDLDALGRCSLADFRALVDAAVVPACRL

>Fitasa-03 (t02_2025_07_03_11_50_3) | T=0.2 | pLDDT=95.85 | Rosetta=-812.332
SEPELRLETVVIVSRHGVRAPTKPTALQKAITPHEWPSWPVKPGTLTPRGAALIASLGAYWR
ERLVALGLLAATGCPAPGEVQILADVDERTRRTGEAFAAGLAPDCDIPVATKPDTTKPDPL
FNPLKTGVCTLDVKAVTDAILAAAGGSMDTFVAARRAEFDLLERVLDFQDSPACLESACRL
TDLYPTELVVTPNNVSLTGAVSLASMLSEIFLLQQAQGLPDPGWGRIRTAEEWDALLALHN
AQFYLLQRTPEVARPRASPLLELIREALTPHPPRPYRYGLTLPTKVLHISGHDTNLANLAG
ALELDWELPGQPDRTPPGGELVFERWRRLADNAWWIQVSAVFQTLQQMRDQTPLSLAAPPA
RVPLTLAGCTDRNARGDCSLADFDALVDAALVPACLL
```

> *Para alineamiento WT vs variantes y conteo exacto de mutaciones, correr el alineamiento pairwise con Biopython o ClustalO usando `1dkl.fa` como query y este archivo como target.*

---

## Anexo B — Inventario de datos

| Directorio (relativo a `~/Escritorio/Protein_Engineering/`) | Contenido clave                                  |
|--------------------------------------------------------------|-------------------------------------------------|
| `1dkl.pdb`                                                   | Estructura base                                 |
| `CONSURF/`                                                   | Resultados completos de ConSurf                 |
| `FIXED_RESIDUES/`                                            | Listas de sitio activo + conservados + set final |
| `MSA-CONSENSUS/`                                             | MSA, consenso, visualización (Conservados.pse)  |
| `HAPYs-Comparison/`                                          | Comparación estructural de HAPs                 |
| `ProteinMPNN-outputs_BN/`                                    | **Corrida final** MPNN + AF3 + top-10 + finales |
| `PoteinMPNN-outputs_FLLS/`                                   | Corrida descartada (variantes de residuos distintos) |
| `Fitasa_PyRosetta/`                                          | Threading, FastRelax, scoring Rosetta, ranking final |

---

*Fin del reporte.*
