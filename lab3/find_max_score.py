import re

# Files to process
files = [
    "bacteria_fungi_blast.txt",
    "bacteria_archaea_blast.txt",
    "fungi_archaea_blast.txt",
    "viral_archaea_blast.txt",
    "viral_bacteria_blast.txt",
    "viral_fungi_blast.txt",
]

# Regular expression to extract the score values
pattern = re.compile(r"Score = ([0-9]+(?:\.[0-9]+)?)\s+bits")

for filename in files:
    max_score = float('-inf')
    max_line = 0
    lineNr = 0

    # Read the file and extract scores
    with open(filename, "r") as file:
        for line in file:
            match = pattern.search(line)
            if match:
                score = float(match.group(1))
                if score > max_score:
                    max_score = score
                    max_line = lineNr
            lineNr += 1

    print("=" * 60)
    print(f"File: {filename}")
    print(f"Maximum Score: {max_score}")
    print("-" * 60)

    # Print context around the max score line
    with open(filename, "r") as file:
        allLines = file.readlines()

    for l in allLines[max_line:max_line + 30]:
        print(l, end="")

    print()  # spacing between files
