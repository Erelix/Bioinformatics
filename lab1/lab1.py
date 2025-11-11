import os
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import product
from collections import Counter
import math
import pandas as pd

def find_valid_orfs(sequence):
    """
    1.  Pateiktoje sekoje fasta formatu surastu visas start ir stop kodonų poras, 
        tarp kurių nebutu stop kodono (ir tiesioginei sekai ir jos reverse komplementui). 
    2.  Kiekvienam stop kodonui parinkti toliausiai nuo jo esanti start kodoną 
        (su salyga, kad tarp ju nera kito stop kodono)
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
    """
    3.  Atfiltruokite visus fragmentus ("tai butu baltymų 
        koduojancios sekos"), kurie trumpesni nei 100 bp.
    """
    return [orf for orf in orfs if (orf[1] - orf[0]) >= min_length_bp]


def translate_orf(sequence, start, stop):
    """
    4.  Konvertuokite koduojancias sekas (start stop kodonu poras) 
        i baltymo seka. Kodonus ir dikodonus analizuokite ne DNR o baltymo 
        lygmenyje (vienas dikodonas - aminorugstis)., t.y tolesnuose 
        zingsniuose - kodonas - viena aminorugstis, dikodonas - dvi.
    """
    dna_seq = Seq(sequence[start:stop])
    prot = str(dna_seq.translate(table=1))
    # Remove trailing stop symbol '*'
    if prot.endswith("*"):
        prot = prot[:-1]
    return prot



def codon_frequencies(sequence):
    """
    5.1 Parasykite funkcijas, kurios ivertintu kodonu daznius 
        (visi imanomi kodonai ir jų atitinkamas daznis  - gali buti 
        nemazai nuliu, jei ju sekoje nerasite).
    """
    sequence = sequence.upper()
    codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3) if len(sequence[i:i+3]) == 3]
    counts = Counter(codons)

    all_codons = [''.join(p) for p in product('ATGC', repeat=3)]
    freqs = {codon: counts.get(codon, 0) for codon in all_codons}
    return freqs


def dicodon_frequencies(sequence):
    """
    5.2 Parasykite funkcijas, kurios ivertintu dikodonu daznius 
        (visi imanomi dikodonai ir jų atitinkamas daznis  - gali buti 
        nemazai nuliu, jei ju sekoje nerasite).
    """
    sequence = sequence.upper()
    codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3) if len(sequence[i:i+3]) == 3]
    dicodons = [codons[i] + codons[i+1] for i in range(len(codons)-1)]
    counts = Counter(dicodons)

    all_dicodons = [''.join(p1 + p2) for p1 in product('ATGC', repeat=3) for p2 in product('ATGC', repeat=3)]
    freqs = {dicodon: counts.get(dicodon, 0) for dicodon in all_dicodons}
    return freqs



def analyze_sequence(seq_record, min_orf_len_bp=100):
    """
    6.  Palyginkite kodonu bei dikodonu daznius tarp visu seku 
    """
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

    # Examples of translated protein sequences (first few)
    # if forward_proteins:
    #     print("\nExample translated protein sequences (forward):")
    #     for i, prot in enumerate(forward_proteins[:3], start=1):
    #         print(f"  ORF{i}: {prot[:60]}{'...' if len(prot) > 60 else ''}")

    # if reverse_proteins:
    #     print("\nExample translated protein sequences (reverse):")
    #     for i, prot in enumerate(reverse_proteins[:3], start=1):
    #         print(f"  ORF{i}: {prot[:60]}{'...' if len(prot) > 60 else ''}")

    # codon_freqs = codon_frequencies(sequence)
    # dicodon_freqs = dicodon_frequencies(sequence)

    # Examples of codon and dicodon frequencies (first few)
    # print(f"Example codon frequencies (first 5): {dict(list(codon_freqs.items())[:5])}")
    # print(f"Example dicodon frequencies (first 5): {dict(list(dicodon_freqs.items())[:5])}")

    print("=" * 50 + "\n")

    return {
        "name": seq_record.id,
        "forward_orfs": forward_orfs,
        "forward_proteins": forward_proteins,
        "reverse_orfs": reverse_orfs,
        "reverse_proteins": reverse_proteins,
    }


# 6. 
def normalize_freqs(freq_dict):
    total = sum(freq_dict.values())
    if total == 0:
        return {k: 0 for k in freq_dict}
    return {k: v / total for k, v in freq_dict.items()}

# 6. √[Σ(freq1_i - freq2_i)²]
def euclidean_distance(vec1, vec2):
    return math.sqrt(sum((vec1[k] - vec2[k]) ** 2 for k in vec1.keys()))


def build_distance_matrix(freqs_dict_all):
    """
    6.  (atstumu matrica - kokia formule naudosite/kaip apskaiciuosite - 
        parasykite ataskaitoje).
    7.  Ivertinkite, ar bakteriniai ir zinduoliu virusai sudaro atskirus 
        klasterius vertinant kodonu/dikodonu dažniu aspektu.
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

# 7.
def simplify_names(df):
    """
    bacterial1 -> B1, mamalian1 -> M1...
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
            new_index.append(idx)

    df.index = new_index
    df.columns = new_index
    return df

# 7.
def save_phylip(matrix, filename):
    n = len(matrix)
    with open(filename, 'w') as f:
        f.write(f"{n}\n")
        for i, row_name in enumerate(matrix.index):
            # PHYLIP: name padded/truncated to 10 characters
            name_fmt = f"{row_name:<10}"[:10]
            row_values = ' '.join(f"{val:.4f}" for val in matrix.iloc[i])
            f.write(f"{name_fmt}{row_values}\n")


def main():
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

            codon_freqs = normalize_freqs(codon_frequencies(seq_str))
            dicodon_freqs = normalize_freqs(dicodon_frequencies(seq_str))

            all_codon_freqs[seq_name] = codon_freqs
            all_dicodon_freqs[seq_name] = dicodon_freqs

            analyze_sequence(seq_record)  # keeps previous analysis

    codon_matrix = build_distance_matrix(all_codon_freqs)
    dicodon_matrix = build_distance_matrix(all_dicodon_freqs)

    codon_matrix = simplify_names(codon_matrix)
    dicodon_matrix = simplify_names(dicodon_matrix)

    codon_matrix = codon_matrix.round(4)
    dicodon_matrix = dicodon_matrix.round(4)

    save_phylip(codon_matrix, "codon_distance_matrix.phy")
    save_phylip(dicodon_matrix, "dicodon_distance_matrix.phy")

    print("\nPHYLIP distance matrices saved as:")
    print("  codon_distance_matrix.phy")
    print("  dicodon_distance_matrix.phy")

if __name__ == "__main__":
    main()
