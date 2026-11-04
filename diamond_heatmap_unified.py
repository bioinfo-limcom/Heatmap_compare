"""
Script to analyze genomic similarities between FASTA files using DIAMOND (6-frame translation).
Can perform all-vs-all analysis in a single folder OR comparative analysis between two folders.
Generates a similarity heatmap and a Neighbor-Joining phylogenetic tree.

================================================================================
DEPENDENCIES & INSTALLATION:
================================================================================
1. Python 3.7+
2. Required Python packages:
   pip install biopython pandas seaborn matplotlib

3. DIAMOND aligner:
   - Download from: https://github.com/bbuchfink/diamond/releases
   - Extract the downloaded archive.
   - Add the executable to your system PATH, OR pass its path via --diamond_path.

================================================================================
USAGE EXAMPLES:
================================================================================
# 1. Single folder (All-vs-All):
python diamond_heatmap_unified.py --folder_a ./genomes --output final_heatmap.png

# 2. Two folders (Folder A vs Folder B):
python diamond_heatmap_unified.py --folder_a ./genomes_A --folder_b ./genomes_B --output final_heatmap.png --diamond_path "C:/path/to/diamond.exe"

# 3. Force re-run (ignoring cached matrices):
python diamond_heatmap_unified.py --folder_a ./genomes --output final_heatmap.png --rerun
"""

import os
import sys
import subprocess
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import shutil
import tempfile
import time
from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

def check_dependencies(diamond_path=None):
    """Check if DIAMOND is available."""
    diamond = "diamond"
    if diamond_path:
        if os.path.isdir(diamond_path):
             diamond = os.path.join(diamond_path, "diamond")
        else:
             diamond = diamond_path
        if not os.path.exists(diamond) and os.path.exists(diamond + ".exe"):
             diamond += ".exe"
    
    if not shutil.which(diamond) and not (os.path.exists(diamond)):
         print(f"Error: 'diamond' not found. Ensure it is installed and in PATH or provided via --diamond_path.")
         print(f"Download it from: https://github.com/bbuchfink/diamond/releases")
         return None
    return diamond

def translate_6_frames(sequence):
    """Translates a DNA sequence into all 6 frames."""
    from Bio.Seq import Seq
    dna = Seq(str(sequence))
    proteins = []
    # Forward frames
    for i in range(3):
        prot = dna[i:].translate(to_stop=False)
        proteins.append((f"f{i+1}", prot))
    # Reverse frames
    rev_dna = dna.reverse_complement()
    for i in range(3):
        prot = rev_dna[i:].translate(to_stop=False)
        proteins.append((f"r{i+1}", prot))
    return proteins

def fragment_sequence(sequence, window, overlap):
    """Generator to yield fragmented sequences."""
    step = window - overlap
    seq_len = len(sequence)
    for i in range(0, seq_len, step):
        yield sequence[i:i+window]
        if i + window >= seq_len: break

def prepare_fragments(input_path, output_path, window=1500, overlap=300):
    """Fragments nucleotide sequences and translates to 6 protein frames."""
    records = []
    for record in SeqIO.parse(input_path, "fasta"):
        seq_len = len(record.seq)
        frags = [record.seq] if seq_len <= window else list(fragment_sequence(record.seq, window, overlap))
        for count, frag_seq in enumerate(frags, 1):
            frag_id = f"{record.id}_frag{count}"
            frames = translate_6_frames(frag_seq)
            for frame_name, prot_seq in frames:
                full_id = f"{frag_id}_{frame_name}"
                records.append(SeqRecord(prot_seq, id=full_id, description=""))
    with open(output_path, "w") as out_handle:
        SeqIO.write(records, out_handle, "fasta")
    return len(records)

def strip_ext(filename):
    """Strips all file extensions from filename."""
    base = str(os.path.basename(filename))
    known_exts = {'.fasta', '.fna', '.fa', '.faa', '.pai', '.pais', '.txt', '.csv', '.nwk', '.png', '.gbk', '.embl', '.gff'}
    
    while True:
        name, ext = os.path.splitext(base)
        ext_str = str(ext).lower()
        if ext_str in known_exts or (ext and len(ext_str) <= 6 and not ext_str.lstrip('.').isdigit()):
            base = name
        else:
            break
    return base

def run_analysis(folder_a, folder_b, output_img, diamond_path=None, rerun=False):
    diamond_bin = check_dependencies(diamond_path)
    if not diamond_bin: sys.exit(1)

    is_two_folders = bool(folder_b)

    exts = ["*.fasta", "*.fa", "*.fna", "*.faa"]
    files_a = sorted([f for ext in exts for f in glob.glob(os.path.join(folder_a, ext))])
    
    if not files_a:
        print(f"Error: No fasta files found in Folder A ({folder_a}).")
        sys.exit(1)

    if is_two_folders:
        files_b = sorted([f for ext in exts for f in glob.glob(os.path.join(folder_b, ext))])
        if not files_b:
            print(f"Error: No fasta files found in Folder B ({folder_b}).")
            sys.exit(1)
        print(f"Mode: Comparative Analysis (Folder A vs Folder B)")
        print(f"Found {len(files_a)} files in Folder A and {len(files_b)} files in Folder B.")
    else:
        files_b = files_a
        print(f"Mode: All-vs-All Analysis (Single Folder)")
        print(f"Found {len(files_a)} files in input folder.")

    print(f"Alignment Method: DIAMOND (6-frame translation)")

    names_a = [strip_ext(f) for f in files_a] # columns
    names_b = [strip_ext(f) for f in files_b] # rows

    base_output = os.path.splitext(output_img)[0]
    cache_path = base_output + '_matrix.csv'
    
    if os.path.exists(cache_path) and not rerun:
        print(f"[DEBUG] Loading cached matrix from {cache_path}")
        matrix = pd.read_csv(cache_path, index_col=0)
    else:
        # Matrix: Rows = B, Columns = A
        matrix = pd.DataFrame(index=names_b, columns=names_a, dtype=float)
        start_time = time.time()
    
        with tempfile.TemporaryDirectory() as temp_dir:
            frag_dir_a = os.path.join(temp_dir, "frag_a")
            os.makedirs(frag_dir_a)
            
            if is_two_folders:
                frag_dir_b = os.path.join(temp_dir, "frag_b")
                os.makedirs(frag_dir_b)
            else:
                frag_dir_b = frag_dir_a

            print(f"\nFragmenting and translating sequences...")
            
            # Prepare Folder A
            frag_map_a = {}
            frag_counts_a = {}
            for f in files_a:
                basename = os.path.basename(f)
                out = os.path.join(frag_dir_a, basename + ".faa")
                frag_map_a[basename] = out
                n = prepare_fragments(f, out)
                frag_counts_a[strip_ext(basename)] = n // 6
                print(f"  [DEBUG] (A) {basename}: {n} protein fragments ({frag_counts_a[strip_ext(basename)]} original)")

            # Prepare Folder B (if different)
            if is_two_folders:
                frag_map_b = {}
                frag_counts_b = {}
                for f in files_b:
                    basename = os.path.basename(f)
                    out = os.path.join(frag_dir_b, basename + ".faa")
                    frag_map_b[basename] = out
                    n = prepare_fragments(f, out)
                    frag_counts_b[strip_ext(basename)] = n // 6
                    print(f"  [DEBUG] (B) {basename}: {n} protein fragments ({frag_counts_b[strip_ext(basename)]} original)")
            else:
                frag_map_b = frag_map_a
                frag_counts_b = frag_counts_a

            print(f"\nRunning DIAMOND alignments...")
            db_dir = os.path.join(temp_dir, "dbs"); os.makedirs(db_dir)

            for name_b_file, path_b in frag_map_b.items():
                db_path = os.path.join(db_dir, name_b_file + ".dmnd")
                subprocess.run([diamond_bin, "makedb", "--in", path_b, "--db", db_path, "--quiet"], check=True)
                
                clean_name_b = strip_ext(name_b_file)
                for name_a_file, path_a in frag_map_a.items():
                    clean_name_a = strip_ext(name_a_file)
                    
                    if not is_two_folders and clean_name_a == clean_name_b:
                        matrix.loc[clean_name_b, clean_name_a] = 100.0
                        continue

                    print(f"  [DEBUG] Aligning {clean_name_a} vs {clean_name_b}...")
                    
                    cmd = [diamond_bin, "blastp", "--query", path_a, "--db", db_path, "--outfmt", "6", "qseqid", "pident", "--max-target-seqs", "1", "--quiet"]
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    
                    if result.returncode != 0 or not result.stdout.strip():
                        matrix.loc[clean_name_b, clean_name_a] = 0.0
                        continue

                    lines = result.stdout.strip().split("\n")
                    best_hits = {}
                    for line in lines:
                        parts = line.split()
                        if len(parts) >= 2:
                            qid, pident = parts[0], float(parts[1])
                            group_id = qid.rsplit('_', 1)[0] # Strip frame
                            if group_id not in best_hits or pident > best_hits[group_id]:
                                best_hits[group_id] = pident
                    
                    total_frags = frag_counts_a[clean_name_a]
                    score = (sum(best_hits.values()) / total_frags) if total_frags > 0 else 0.0
                    matrix.loc[clean_name_b, clean_name_a] = score
                    print(f"    [DEBUG] Result: {score:.1f}% ({len(best_hits)}/{total_frags} fragments hit)")

        end_time = time.time()
        print(f"\nTotal analysis time: {end_time - start_time:.2f} seconds")
        matrix = matrix.apply(pd.to_numeric).fillna(0.0)
        matrix.to_csv(cache_path)

    # ---- Build Tree and Heatmap ----
    def build_nj_tree(item_names, data_matrix, is_rows=True):
        if len(item_names) < 3:
            return None, sorted(item_names), {}
        
        dist_data = []
        for i, s1 in enumerate(item_names):
            row_dists = []
            for j, s2 in enumerate(item_names):
                if i == j: row_dists.append(0.0)
                elif j < i:
                    v1 = data_matrix.loc[s1].values if is_rows else data_matrix[s1].values
                    v2 = data_matrix.loc[s2].values if is_rows else data_matrix[s2].values
                    d = sum(abs(v1 - v2)) / len(v1)
                    row_dists.append(d)
            dist_data.append(row_dists)
        
        clean_names = [s.replace('.', '_').replace('-', '_').replace(' ', '_') for s in item_names]
        name_map = dict(zip(clean_names, item_names))
        
        constructor = DistanceTreeConstructor()
        matrix_dist = DistanceMatrix(clean_names, dist_data)
        tree = constructor.nj(matrix_dist)
        
        avg_profile_val = data_matrix.mean(axis=1) if is_rows else data_matrix.mean(axis=0)
        outgroup_raw = avg_profile_val.idxmin()
        outgroup_clean = outgroup_raw.replace('.', '_').replace('-', '_').replace(' ', '_')
        outgroup_clades = [c for c in tree.get_terminals() if c.name == outgroup_clean]
        if outgroup_clades: tree.root_with_outgroup(outgroup_clades[0])
        
        leaf_order = [name_map.get(c.name, c.name) for c in tree.get_terminals()]
        return tree, leaf_order, name_map

    print("\n[DEBUG] Building trees for the matrix...")
    
    # Tree for rows (Folder B / Main Genomes)
    row_tree, row_order, row_map = build_nj_tree(list(matrix.index), matrix, is_rows=True)
    
    # Tree for columns (Folder A)
    if not is_two_folders or set(matrix.index) == set(matrix.columns):
        col_tree, col_order, col_map = None, row_order, row_map
    else:
        col_tree, col_order, col_map = build_nj_tree(list(matrix.columns), matrix, is_rows=False)

    if row_tree:
        Phylo.write(row_tree, base_output + '_row_tree.nwk', 'newick')
    if col_tree:
        Phylo.write(col_tree, base_output + '_col_tree.nwk', 'newick')
    
    heatmap_out = base_output + "_heatmap.png"
    tree_out = base_output + "_phylogenetic_tree.png"
    
    if row_tree:
        print(f"\n[DEBUG] Generating Phylogenetic Tree (Rows): {tree_out}")
        tree_height = max(len(row_order) * 0.4, 6)
        fig_tree, ax_tree = plt.subplots(figsize=(8, tree_height))
        Phylo.draw(row_tree, axes=ax_tree, do_show=False, label_func=lambda c: row_map.get(c.name, c.name) if c.is_terminal() else '')
        
        ax_tree.set_title('Phylogenetic Tree', fontsize=14, fontweight='bold', pad=20)
        ax_tree.spines['top'].set_visible(False)
        ax_tree.spines['right'].set_visible(False)
        ax_tree.spines['bottom'].set_visible(True)
        ax_tree.set_xlabel('Substitution per fragment', fontsize=10, style='italic')
        plt.tight_layout()
        plt.savefig(tree_out, dpi=300, bbox_inches='tight')
        plt.close()

    if col_tree:
        col_tree_out = base_output + "_col_tree.png"
        print(f"[DEBUG] Generating Phylogenetic Tree (Cols): {col_tree_out}")
        tree_h = max(len(col_order) * 0.4, 6)
        fig_c, ax_c = plt.subplots(figsize=(8, tree_h))
        Phylo.draw(col_tree, axes=ax_c, do_show=False, label_func=lambda c: col_map.get(c.name, c.name) if c.is_terminal() else '')
        ax_c.set_title('Column Phylogenetic Tree', fontsize=14, fontweight='bold', pad=20)
        ax_c.spines['top'].set_visible(False)
        ax_c.spines['right'].set_visible(False)
        ax_c.set_xlabel('Substitution per fragment', fontsize=10, style='italic')
        plt.tight_layout()
        plt.savefig(col_tree_out, dpi=300, bbox_inches='tight')
        plt.close()

    print(f"[DEBUG] Generating Heatmap: {heatmap_out}")
    matrix_ordered = matrix.loc[row_order][col_order]
    
    cell_size = 0.6
    fig_w = max(len(matrix.columns) * cell_size + 2, 10)
    fig_h = max(len(matrix.index) * cell_size + 2, 8)
    
    plt.figure(figsize=(fig_w, fig_h))
    
    font_size = int(cell_size * 72 * 0.3) 
    font_size = min(max(font_size, 6), 14)

    ax_heat = sns.heatmap(matrix_ordered, annot=True, 
                          cmap=sns.blend_palette(["green", "yellow", "orange", "red"], as_cmap=True), 
                          fmt=".0f", linewidths=0.5, linecolor='gray', 
                          annot_kws={"size": font_size, "weight": "bold"})
    
    title = 'Genomic Similarity Matrix (%)'
    if is_two_folders:
        title += f'\n(Rows: Folder B, Cols: Folder A)'
        
    ax_heat.set_title(title, fontsize=16, fontweight='bold', pad=20)
    ax_heat.set_xlabel('', fontsize=12, fontweight='bold')
    ax_heat.set_ylabel('', fontsize=12, fontweight='bold')
    
    plt.savefig(heatmap_out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n[DONE] Finished generating files. Output images saved to:")
    print(f"  - {heatmap_out}")
    if row_tree:
        print(f"  - {tree_out}")
    if col_tree:
        print(f"  - {base_output}_col_tree.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze genomic similarities using DIAMOND (6-frame translation).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single folder (All-vs-All):
  python %(prog)s --folder_a ./genomes --output final_heatmap.png
  
  # Two folders (Folder A vs Folder B):
  python %(prog)s --folder_a ./genomes_A --folder_b ./genomes_B --output final_heatmap.png
        """
    )
    parser.add_argument("-a", "--folder_a", required=True, help="Path to first folder of FASTA files (or single folder for all-vs-all)")
    parser.add_argument("-b", "--folder_b", help="Path to second folder of FASTA files (optional)", default=None)
    parser.add_argument("-o", "--output", required=True, help="Output path for heatmap image (e.g. final_heatmap.png)")
    parser.add_argument("--diamond_path", help="Path to DIAMOND executable", default=None)
    parser.add_argument("--rerun", action="store_true", help="Force re-run (ignore cached matrices)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.folder_a):
        print(f"Error: Folder A '{args.folder_a}' does not exist.")
        sys.exit(1)
        
    if args.folder_b and not os.path.exists(args.folder_b):
        print(f"Error: Folder B '{args.folder_b}' does not exist.")
        sys.exit(1)
        
    run_analysis(args.folder_a, args.folder_b, args.output, args.diamond_path, args.rerun)
