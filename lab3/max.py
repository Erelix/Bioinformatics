import re

# Files to process
files = [
    "bacteria_fungi_blast.txt",
    "bacteria_archaea_blast.txt",
    "fungi_archaea_blast.txt",
]

# Patterns
score_pattern = re.compile(r"Score = ([0-9]+(?:\.[0-9]+)?)\s+bits")
sbjct_pattern = re.compile(r"Sbjct\s+\d+\s+([A-Z\-]+)\s+\d+")

# Output FASTA file
with open("max_sequences.fasta", "w") as out_file:

    for filename in files:
        with open(filename) as f:
            lines = f.readlines()

        # Find max score
        max_score = 0
        max_line_nr = 0
        for i, line in enumerate(lines):
            m = score_pattern.search(line)
            if m:
                score = float(m.group(1))
                if score > max_score:
                    max_score = score
                    max_line_nr = i

        # Extract full Sbjct sequence for max hit
        seq_lines = []
        i = max_line_nr
        while i < len(lines):
            line = lines[i]
            if line.startswith("Sbjct"):
                m = sbjct_pattern.search(line)
                if m:
                    seq_lines.append(m.group(1).replace("-", ""))
            elif line.startswith("Score") and i != max_line_nr:
                # Stop at next hit
                break
            i += 1

        full_seq = "".join(seq_lines)

        # Write to FASTA
        out_file.write(f">{filename}|Score={max_score}\n")
        # Wrap sequence at 80 characters per line
        for j in range(0, len(full_seq), 80):
            out_file.write(full_seq[j:j+80] + "\n")
