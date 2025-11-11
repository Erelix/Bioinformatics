import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import product
from collections import Counter
import math
import pandas as pd

def find_valid_orfs(sequence):
    """
    For each reading frame, find stop codons and for each stop codon choose
    the farthest start codon (ATG) before it, provided there is no other stop
    codon in between.
    Returns list of tuples (start, stop, frame).
    """
    valid_orfs = []
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}

    for frame in range(3):
        frame_seq = sequence[frame:]
        codons = [frame_seq[i:i+3] for i in range(0, len(frame_seq)-2, 3)]
        start_positions = [i*3 + frame for i, c in enumerate(codons) if c == start_codon]
        stop_positions = [i*3 + frame for i, c in enumerate(codons) if c in stop_codons]

        for stop_idx in stop_positions:
            farthest_start = None
            for start_idx in reversed(start_positions):
                if start_idx < stop_idx:
                    if any(s > start_idx and s < stop_idx for s in stop_positions):
                        break
                    farthest_start = start_idx
                    break
            if farthest_start is not None:
                valid_orfs.append((farthest_start, stop_idx + 3, frame))
    return valid_orfs


def filter_orfs_min_length(orfs, min_length_bp=100):
    """Keep only ORFs >= min_length_bp."""
    return [orf for orf in orfs if (orf[1] - orf[0]) >= min_length_bp]


def translate_orf(sequence, start, stop):
    """Translate DNA to protein using Biopython codon table 1."""
    dna_seq = Seq(sequence[start:stop])
    prot = str(dna_seq.translate(table=1))
    # Remove trailing stop symbol '*'
    if prot.endswith("*"):
        prot = prot[:-1]
    return prot



def codon_frequencies(sequence):
    """
    Returns a dictionary of codon frequencies for all 64 codons (AAA..TTT),
    even if they don't appear in the sequence.
    """
    sequence = sequence.upper()
    codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3) if len(sequence[i:i+3]) == 3]
    counts = Counter(codons)

    all_codons = [''.join(p) for p in product('ATGC', repeat=3)]
    freqs = {codon: counts.get(codon, 0) for codon in all_codons}
    return freqs


def dicodon_frequencies(sequence):
    """
    Returns a dictionary of dicodon (6-base) frequencies for all 4096 possible
    dicodons (AAA AAA .. TTT TTT), even if they don't appear.
    """
    sequence = sequence.upper()
    codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3) if len(sequence[i:i+3]) == 3]
    dicodons = [codons[i] + codons[i+1] for i in range(len(codons)-1)]
    counts = Counter(dicodons)

    all_dicodons = [''.join(p1 + p2) for p1 in product('ATGC', repeat=3) for p2 in product('ATGC', repeat=3)]
    freqs = {dicodon: counts.get(dicodon, 0) for dicodon in all_dicodons}
    return freqs



def analyze_sequence(seq_record, min_orf_len_bp=100):
    """Find valid ORFs, filter short ones, and translate to protein."""
    sequence = str(seq_record.seq).upper()
    print(f"Analyzing: {seq_record.id}")
    print(f"Sequence length: {len(sequence)}")

    # Forward strand
    forward_orfs = filter_orfs_min_length(find_valid_orfs(sequence), min_orf_len_bp)
    forward_proteins = [translate_orf(sequence, s, e) for s, e, f in forward_orfs]
    print(f"Forward ORFs (≥{min_orf_len_bp}bp): {len(forward_orfs)}")

    # Reverse strand
    rev_comp = str(seq_record.seq.reverse_complement()).upper()
    reverse_orfs = filter_orfs_min_length(find_valid_orfs(rev_comp), min_orf_len_bp)
    reverse_proteins = [translate_orf(rev_comp, s, e) for s, e, f in reverse_orfs]
    print(f"Reverse ORFs (≥{min_orf_len_bp}bp): {len(reverse_orfs)}")

    # Show sample translated protein sequences (first few)
    if forward_proteins:
        print("\nExample translated protein sequences (forward):")
        for i, prot in enumerate(forward_proteins[:3], start=1):
            print(f"  ORF{i}: {prot[:60]}{'...' if len(prot) > 60 else ''}")

    if reverse_proteins:
        print("\nExample translated protein sequences (reverse):")
        for i, prot in enumerate(reverse_proteins[:3], start=1):
            print(f"  ORF{i}: {prot[:60]}{'...' if len(prot) > 60 else ''}")

    codon_freqs = codon_frequencies(sequence)
    dicodon_freqs = dicodon_frequencies(sequence)

    print(f"Example codon frequencies (first 10): {dict(list(codon_freqs.items())[:10])}")
    print(f"Example dicodon frequencies (first 5): {dict(list(dicodon_freqs.items())[:5])}")


    print("=" * 50 + "\n")

    return {
        "name": seq_record.id,
        "forward_orfs": forward_orfs,
        "forward_proteins": forward_proteins,
        "reverse_orfs": reverse_orfs,
        "reverse_proteins": reverse_proteins,
    }

def normalize_freqs(freq_dict):
    total = sum(freq_dict.values())
    if total == 0:
        return {k: 0 for k in freq_dict}
    return {k: v / total for k, v in freq_dict.items()}


def euclidean_distance(vec1, vec2):
    return math.sqrt(sum((vec1[k] - vec2[k]) ** 2 for k in vec1.keys()))


def build_distance_matrix(freqs_dict_all):
    """
    freqs_dict_all: dict {sequence_name: {codon: frequency}}
    returns: pandas DataFrame (distance matrix)
    """
    seq_names = list(freqs_dict_all.keys())
    matrix = []

    for name1 in seq_names:
        row = []
        for name2 in seq_names:
            dist = euclidean_distance(freqs_dict_all[name1], freqs_dict_all[name2])
            row.append(dist)
        matrix.append(row)

    return pd.DataFrame(matrix, index=seq_names, columns=seq_names)    


def simplify_names(df):
    """
    Rename long sequence names to short labels:
    bacterial1 -> B1, mamalian1 -> M1, etc.
    """
    rename_map = {
        'bacterial1': 'B1',
        'bacterial2': 'B2',
        'bacterial3': 'B3',
        'bacterial4': 'B4',
        'mamalian1': 'M1',
        'mamalian2': 'M2',
        'mamalian3': 'M3',
        'mamalian4': 'M4',
    }

    new_index = []
    for idx in df.index:
        for old, new in rename_map.items():
            if old in idx:
                new_index.append(new)
                break
        else:
            new_index.append(idx)  # fallback if no match

    df.index = new_index
    df.columns = new_index
    return df


def main_multiple_files():
    fasta_files = {
        'bacterial1': r".\viruses\viruses\data\bacterial1.fasta",
        'bacterial2': r".\viruses\viruses\data\bacterial2.fasta",
        'bacterial3': r".\viruses\viruses\data\bacterial3.fasta",
        'bacterial4': r".\viruses\viruses\data\bacterial4.fasta",
        'mamalian1': r".\viruses\viruses\data\mamalian1.fasta",
        'mamalian2': r".\viruses\viruses\data\mamalian2.fasta",
        'mamalian3': r".\viruses\viruses\data\mamalian3.fasta",
        'mamalian4': r".\viruses\viruses\data\mamalian4.fasta",
    }

    all_codon_freqs = {}
    all_dicodon_freqs = {}

    for file_name, file_path in fasta_files.items():
        if not os.path.exists(file_path):
            print(f"Warning: File not found - {file_path}")
            continue

        for seq_record in SeqIO.parse(file_path, "fasta"):
            seq_name = f"{file_name}_{seq_record.id}"
            seq_str = str(seq_record.seq).upper()

            # Codon & dicodon frequencies (normalized)
            codon_freqs = normalize_freqs(codon_frequencies(seq_str))
            dicodon_freqs = normalize_freqs(dicodon_frequencies(seq_str))

            all_codon_freqs[seq_name] = codon_freqs
            all_dicodon_freqs[seq_name] = dicodon_freqs

            analyze_sequence(seq_record)  # keep your previous analysis

    # Build distance matrices
    codon_matrix = build_distance_matrix(all_codon_freqs)
    dicodon_matrix = build_distance_matrix(all_dicodon_freqs)

    # Simplify names (B1, M1, etc.)
    codon_matrix = simplify_names(codon_matrix)
    dicodon_matrix = simplify_names(dicodon_matrix)

    # Round to 4 decimal places
    codon_matrix = codon_matrix.round(4)
    dicodon_matrix = dicodon_matrix.round(4)

    print("\n=== CODON DISTANCE MATRIX ===")
    print(codon_matrix)

    print("\n=== DICODON DISTANCE MATRIX ===")
    print(dicodon_matrix)

    with open("codon_distance_matrix.txt", "w") as f:
        f.write("=== CODON DISTANCE MATRIX ===\n")
        f.write(codon_matrix.to_string(float_format="{:.4f}".format))
        f.write("\n")

    with open("dicodon_distance_matrix.txt", "w") as f:
        f.write("=== DICODON DISTANCE MATRIX ===\n")
        f.write(dicodon_matrix.to_string(float_format="{:.4f}".format))
        f.write("\n")




if __name__ == "__main__":
    main_multiple_files()
