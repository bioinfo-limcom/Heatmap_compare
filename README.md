# 🧬 Heatmap_compare

A fast and versatile Python pipeline to analyze genomic similarities between FASTA files using **DIAMOND** (6-frame translation).  
Supports **all-vs-all** comparisons within a single directory or **comparative** analyses between two distinct sets of genomes/regions (Folder A vs. Folder B), automatically generating a similarity heatmap and Neighbor-Joining (NJ) phylogenetic trees.

---

## 📌 Overview

Assessing sequence similarity and conservation across bacterial genomes, genomic islands, plasmids, or mobile genetic elements often requires intensive all-vs-all alignments.  
**Heatmap_compare** (`diamond_heatmap_unified.py`) streamlines this process by:
1. Fragmenting DNA sequences into overlapping windows (1,500 bp with 300 bp overlap).
2. Translating nucleotide fragments in all 6 reading frames.
3. Conducting ultra-fast protein alignments using **DIAMOND**.
4. Computing a normalized similarity percentage matrix.
5. Building and plotting Neighbor-Joining (NJ) phylogenetic trees.
6. Generating publication-ready annotated heatmaps (300 DPI).

---

## 🛠️ Features

- **Dual Comparison Modes**:
  - **All-vs-All** within a single folder.
  - **Comparative Analysis** between two separate directories (Folder A vs. Folder B).
- **6-Frame Translation**: Translates DNA fragments across all six reading frames for robust protein-level alignment.
- **Ultra-Fast Alignments**: Leverages the high performance of DIAMOND.
- **Phylogenetic Reconstruction**: Automatically builds Neighbor-Joining (NJ) distance trees rooted via an automatic outgroup algorithm.
- **Newick Tree Export**: Saves phylogenetic trees in Newick format (`.nwk`) for easy inspection in iTOL, FigTree, or MEGA.
- **Publication-Ready Figures**: Outputs high-resolution (300 DPI) heatmaps with customizable cell scaling and clean phylogenetic tree diagrams.
- **Smart Caching**: Saves computed distance matrices as `.csv`. Re-runs with different plotting parameters load instantly unless `--rerun` is explicitly passed.

---

## 📦 Requirements & Dependencies

### 1. Python Environment
- Python 3.7+
- Required libraries:
  - `biopython`
  - `pandas`
  - `seaborn`
  - `matplotlib`

### 2. External Tools
- **DIAMOND Aligner** (v2.0+ recommended)
  - Download from: [https://github.com/bbuchfink/diamond/releases](https://github.com/bbuchfink/diamond/releases)
  - Add `diamond` to your system `PATH`, or provide the binary path via `--diamond_path`.

---

## 🚀 Installation

Clone the repository and install the Python dependencies:

```bash
git clone https://github.com/bioinfo-limcom/Heatmap_compare.git
cd Heatmap_compare

pip install biopython pandas seaborn matplotlib
```

Ensure DIAMOND is installed and accessible:
```bash
diamond --version
```

---

## ▶️ Usage

### 1. Single Folder (All-vs-All Analysis)
Compare all FASTA files in a single directory against each other:

```bash
python diamond_heatmap_unified.py \
    --folder_a ./genomes \
    --output results_heatmap.png
```

### 2. Two Folders (Comparative Analysis: Folder A vs. Folder B)
Compare files in Folder A against files in Folder B:

```bash
python diamond_heatmap_unified.py \
    --folder_a ./genomic_islands \
    --folder_b ./reference_genomes \
    --output comparison_heatmap.png
```

### 3. Specify DIAMOND Executable Path
If DIAMOND is not in your system `PATH`:

```bash
python diamond_heatmap_unified.py \
    --folder_a ./genomes \
    --output results_heatmap.png \
    --diamond_path /usr/local/bin/diamond
```
*(On Windows: `--diamond_path "C:\path\to\diamond.exe"`)*

### 4. Force Re-run (Bypass Cache)
Recompute alignments from scratch, ignoring existing cached `*_matrix.csv`:

```bash
python diamond_heatmap_unified.py \
    --folder_a ./genomes \
    --output results_heatmap.png \
    --rerun
```

---

## ⚙️ Parameters & Arguments

| Argument | Shorthand | Required | Description |
| :--- | :---: | :---: | :--- |
| `--folder_a` | `-a` | **Yes** | Path to the primary folder containing FASTA files (`.fasta`, `.fa`, `.fna`, `.faa`). |
| `--folder_b` | `-b` | No | Path to the second folder of FASTA files for two-way comparison. |
| `--output` | `-o` | **Yes** | Base output path for the heatmap image (e.g., `heatmap.png`). |
| `--diamond_path` | | No | Explicit path to the DIAMOND executable or its containing directory. |
| `--rerun` | | No | Force recalculation of similarity matrix, bypassing cached `.csv`. |

---

## 🔄 Workflow Pipeline

```
[ Nucleotide FASTA Files ]
           │
           ▼
[ Sliding Window Fragmentation ]  (1500 bp window, 300 bp overlap)
           │
           ▼
[ 6-Frame In Silico Translation ]
           │
           ▼
[ DIAMOND DB & BLASTP Alignment ]
           │
           ▼
[ Normalized Similarity Matrix ] ──► Saved as *_matrix.csv
           │
     ┌─────┴─────────────────────┐
     ▼                           ▼
[ NJ Tree Construction ]    [ Heatmap Generation ]
  - Pairwise Distance         - Color-graded similarity (%)
  - Outgroup rooting          - Dynamic dimension scaling
  - Exported as *.nwk / *.png - Exported as 300 DPI *.png
```

---

## 📊 Output Files

Given an output name like `--output results.png`, the pipeline produces:

- **`results_heatmap.png`**: High-resolution heatmap visualizing pairwise similarity percentages.
- **`results_matrix.csv`**: Tabular CSV matrix containing the raw pairwise similarity scores.
- **`results_phylogenetic_tree.png`**: Neighbor-Joining phylogenetic tree plot (row entities).
- **`results_row_tree.nwk`**: Newick formatted tree file for downstream tree visualization software.
- **`results_col_tree.png`** / **`results_col_tree.nwk`**: Column phylogenetic tree (generated during two-folder comparison when sets differ).

---

## 📁 Expected Directory Structure

Supported file extensions for inputs include `.fasta`, `.fa`, `.fna`, and `.faa`.

### Example (Single Folder Mode):
```
genomes/
├── Strain_01.fasta
├── Strain_02.fasta
└── Strain_03.fasta
```

### Example (Two Folder Mode):
```
dataset/
├── query_islands/
│   ├── Island_A.fasta
│   └── Island_B.fasta
└── reference_genomes/
    ├── Ref_Strain_01.fna
    └── Ref_Strain_02.fna
```

---

## ❗ Notes & Best Practices

- **Sequence Formats**: Nucleotide FASTA files (`.fasta`, `.fa`, `.fna`) are fragmented and translated in 6 frames.
- **Cache Acceleration**: Once an analysis is run, the script creates `*_matrix.csv`. If you rerun with the same output prefix to adjust figure parameters, the alignment step will be skipped automatically unless `--rerun` is supplied.
- **File Names**: The script automatically strips file extensions when labeling heatmap rows and columns for clean presentation.

---

## 📄 License

This project is freely available for academic and scientific use.  
For commercial use, please contact the authors.

---

## ✨ Author

**Andrei Giacchetto Felice**  
Laboratory of Immunology and Omics Sciences (LimCom)  
Universidade Federal do Triângulo Mineiro (UFTM)
