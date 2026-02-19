import re

files = ["viral_archaea_blast.txt", "viral_bacteria_blast.txt", "viral_fungi_blast.txt"]
top_n_per_file = 3  
max_sequences = 6   
output_fasta = "top_hits_for_alphafold.fasta"

def parse_blast(file):
    hits = []
    with open(file) as f:
        content = f.read()
    
    # Split by hits (">" lines)
    for block in content.split(">")[1:]:
        header_line = block.splitlines()[0]
        score_match = re.search(r"Score\s*=\s*([\d.]+)\s*bits", block)
        if score_match:
            score = float(score_match.group(1))
            hits.append((score, header_line, block))
    
    # Sort by score descending
    hits.sort(reverse=True, key=lambda x: x[0])
    return hits[:top_n_per_file]

def extract_sequence(block):
    seq = ""
    for line in block.splitlines():
        if line.startswith("Query"):
            parts = line.split()
            seq += parts[2].replace("-", "")
    return seq

seen_sequences = set()
total_sequences = 0

with open(output_fasta, "w") as out:
    for file in files:
        top_hits = parse_blast(file)
        for score, header, block in top_hits:
            seq = extract_sequence(block)
            if seq not in seen_sequences:
                seen_sequences.add(seq)
                fasta_header = f">{header.strip()} | Score={score} | Source={file}"
                out.write(f"{fasta_header}\n{seq}\n")
                total_sequences += 1
                if total_sequences >= max_sequences:
                    break
        if total_sequences >= max_sequences:
            break

print(f"Saved {total_sequences} unique sequences to {output_fasta}")
