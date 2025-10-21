import re
import os
from Bio import SeqIO
from Bio.Seq import Seq

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

    print("=" * 50 + "\n")

    return {
        "name": seq_record.id,
        "forward_orfs": forward_orfs,
        "forward_proteins": forward_proteins,
        "reverse_orfs": reverse_orfs,
        "reverse_proteins": reverse_proteins,
    }


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

    for file_name, file_path in fasta_files.items():
        if not os.path.exists(file_path):
            print(f"Warning: File not found - {file_path}")
            continue

        print("=" * 60)
        print(f"File: {file_name}.fasta")
        print("=" * 60)

        for seq_record in SeqIO.parse(file_path, "fasta"):
            analyze_sequence(seq_record)


if __name__ == "__main__":
    main_multiple_files()
