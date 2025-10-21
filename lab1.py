import re
import os
from Bio import SeqIO
from Bio.Seq import Seq

def find_valid_orfs(sequence):
    """
    1. Pateiktoje sekoje fasta formatu surastu visas start ir stop kodonų poras, tarp kurių nebutu stop kodono (ir tiesioginei sekai ir jos reverse komplementui).
    2. Kiekvienam stop kodonui parinkti toliausiai nuo jo esanti start kodoną (su salyga, kad tarp ju nera kito stop kodono)
    """
    valid_orfs = []
    
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']

    for frame in range(3):
        frame_seq = sequence[frame:]
        
        # Find all start codons in this frame
        start_positions = []
        for match in re.finditer(start_codon, frame_seq):
            start_positions.append(match.start())
        
        # For each start codon, find the next stop codon
        for start_idx in start_positions:
            absolute_start = frame + start_idx
            remaining_seq = frame_seq[start_idx:]
            
            # Check every codon after the start
            found_stop = False
            
            for i in range(3, len(remaining_seq) - 2, 3):
                codon = remaining_seq[i:i+3]
                
                if codon in stop_codons:
                    # Found a stop codon - valid ORF
                    absolute_stop = absolute_start + i + 3
                    valid_orfs.append((absolute_start, absolute_stop, frame))
                    found_stop = True
                    break
    
    return valid_orfs

def analyze_sequence(seq_record):
    """Analyze a sequence and its reverse complement for valid ORFs"""
    sequence = str(seq_record.seq).upper()
    
    print(f"Analyzing: {seq_record.id}")
    print(f"Sequence length: {len(sequence)}")
    
    forward_orfs = find_valid_orfs(sequence)
    print(f"Forward strand - {len(forward_orfs)} ORFs")
    
    rev_comp = str(seq_record.seq.reverse_complement())
    reverse_orfs = find_valid_orfs(rev_comp)
    print(f"Reverse complement - {len(reverse_orfs)} ORFs")

    print(f"{'='*40}"+"\n")
    
    # Print details of found ORFs (first 5 only)
    # if forward_orfs:
    #     print("\nForward strand ORFs (first 5):")
    #     for start, stop, frame in forward_orfs[:5]:
    #         length = stop - start
    #         orf_seq = sequence[start:stop]
    #         print(f"  Frame {frame}: Start={start}, Stop={stop}, Length={length}")
    #         print(f"    Start: {orf_seq[:10]}... Stop: ...{orf_seq[-10:]}")
    
    # if reverse_orfs:
    #     print("\nReverse complement ORFs (first 5):")
    #     for start, stop, frame in reverse_orfs[:5]:
    #         length = stop - start
    #         orf_seq = rev_comp[start:stop]
    #         print(f"  Frame {frame}: Start={start}, Stop={stop}, Length={length}")
    #         print(f"    Start: {orf_seq[:10]}... Stop: ...{orf_seq[-10:]}")
    
    return {
        'name': seq_record.id,
        'forward_orfs': forward_orfs,
        'reverse_orfs': reverse_orfs,
        'total_orfs': len(forward_orfs) + len(reverse_orfs)
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
    
    print("Checking files:", fasta_files)
    
    all_results = {}
    
    for file_name, file_path in fasta_files.items():
        if not os.path.exists(file_path):
            print(f"Warning: File not found - {file_path}")
            continue
        
        try:
            print(f"{'='*40}")
            print(f"File: {file_name}"+".fasta")
            
            sequences = list(SeqIO.parse(file_path, "fasta"))
            file_results = []
            
            for seq_record in sequences:
                result = analyze_sequence(seq_record)
                file_results.append(result)
            
            all_results[file_name] = file_results
            
        except Exception as e:
            print(f"Error processing {file_name}: {e}")
    
    # # Final summary comparing all files
    # print("\n" + "=" * 60)
    # print("COMPARATIVE SUMMARY")
    # print("=" * 60)
    # for file_name, results in all_results.items():
    #     total_orfs = sum(result['total_orfs'] for result in results)
    #     total_forward = sum(len(result['forward_orfs']) for result in results)
    #     total_reverse = sum(len(result['reverse_orfs']) for result in results)
        
    #     print(f"{file_name}:")
    #     print(f"  Total ORFs: {total_orfs}")
    #     print(f"  Forward: {total_forward}, Reverse: {total_reverse}")
    #     print(f"  Sequences analyzed: {len(results)}")


if __name__ == "__main__":
    main_multiple_files()