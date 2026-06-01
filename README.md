# linkerology

A cheminformatics pipeline for automated linker classification, physicochemical profiling,
conformational analysis, and interactive reporting of PROTAC and bifunctional degrader linkers.
Developed to support the 2026 Linkerology Review in *Trends in Biochemical Sciences*.

---

## Features

### Automated linker classification
Rule-based SMARTS matching for eight recurrent chemotype families:
- Aliphatic (acyclic sp³ chains)
- Ether / polyether (including PEG motifs)
- Amide / amine
- Alicyclic / heterocyclic (piperazine, morpholine, piperidine…)
- Aromatic / heteroaromatic
- Triazole (click chemistry)
- Rigid (alkynes, alkenes)
- Functionalised (disulfides, azo, azide, carbonate…)

Linkers matching more than one motif are automatically assigned a `Hybrid(A + B + …)` label
and are expanded for all downstream analyses.

### Physicochemical property calculation
Heavy atoms, rotatable bonds, TPSA, cLogP, fraction sp³, total rings, aromatic rings.

### Linker length
Estimated as the **maximum** heavy-atom end-to-end distance across an ensemble of 30
energy-minimised conformers (RDKit ETKDGv3 + UFF force field, `randomSeed = 42`).

### Structured output
All outputs are written to the `reports/` directory:
- `classified_linkers.csv`
- `classified_linkers.xlsx` (includes embedded molecular structure images)
- `category_counts.csv`
- `data_by_category.csv` (hybrid-expanded)
- `category_map.csv` (binary presence/absence per chemotype)
- `summary.html` (interactive Plotly report with clickable category filtering)

---

## Installation

```bash
conda create -n linker-env python=3.11 -y
conda activate linker-env
conda install -c conda-forge rdkit pandas numpy openpyxl plotly scikit-learn tqdm pillow matplotlib -y
```

> **Note:** For static SVG export from `4_chord.py`, also install `kaleido`:
> `pip install kaleido`

---

## Usage

### 1. Classify and profile linkers
```bash
python 1_linkeranalysis.py your_linkers.sdf
```
Outputs:
- `classified_linkers.xlsx` — in the working directory, with embedded structure images
- `reports/classified_linkers.csv`, `reports/category_counts.csv`, `reports/data_by_category.csv`,
  `reports/category_map.csv`, `reports/summary.html`

### 2. PCA
```bash
python 2_cluster.py classified_linkers.xlsx
```
Outputs `classified_linkers_PCA.csv` with PC1/PC2 coordinates and chemotype labels,
ready to import into GraphPad Prism or any plotting tool.

### 3. Linker length (end-to-end distance)
```bash
python 3_length.py classified_linkers.xlsx
```
Generates `classified_linkers_distances.xlsx` with `dist_max_A` (maximum heavy-atom
end-to-end distance across 30 ETKDGv3 + UFF conformers, Å) for each linker.

Options:
```
--nconfs INT    Number of conformers per molecule (default: 30)
--workers INT   Parallel worker processes (default: min(cpu_count, 8))
--output PATH   Custom output file path
```

### 4. Chord diagram of chemotype co-occurrence
```bash
python 4_chord.py
# or specify paths explicitly:
python 4_chord.py --input reports/category_map.csv --outdir graphs/
# with optional static SVG export (requires kaleido):
python 4_chord.py --svg
```
Chord width is proportional to log-scaled pairwise co-occurrence frequency.

---

## Recommended workflow

```
1_linkeranalysis.py  →  classified_linkers.xlsx  →  2_cluster.py  →  classified_linkers_PCA.csv
                                                 →  3_length.py   →  classified_linkers_distances.xlsx
                     →  reports/category_map.csv →  4_chord.py    →  graphs/chord_cooccurrence.html
```

---

## Input data

The published analysis uses linker structures from **PROTAC-DB**
(https://cadd.zju.edu.cn/protacdb/, accessed May 2025), operated by Tingjun Hou's
Group at Zhejiang University. PROTAC-DB data is subject to their terms of use;
users must download it directly from the platform and may not redistribute it.

To reproduce:
1. Go to https://cadd.zju.edu.cn/protacdb/ → Linker section → select all entries → export as SDF.
2. Run `python 1_linkeranalysis.py linkers.sdf`.

The pipeline accepts any standard SDF file — the `_Name` property is used as the
linker identifier. Conformer generation in `3_length.py` uses `randomSeed = 42`.

## Notes
- Hybrid categories are automatically expanded for all downstream analyses.
- Designed for datasets of hundreds to thousands of linkers.

---


## Authors
Rubén Prieto-Díaz (2025)
Centre for Targeted Protein Degradation (CeTPD), University of Dundee
CiQUS – University of Santiago de Compostela

---

## License
Released under the MIT License. You are free to use, modify, and distribute this software
with proper attribution.
