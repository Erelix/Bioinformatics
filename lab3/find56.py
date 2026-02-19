import re

# Files to check
blast_files = ["viral_archaea_blast.txt", "viral_bacteria_blast.txt", "viral_fungi_blast.txt"]

top_hits_per_file = 2  # get 2 max from each file

result = []

for bf in blast_files:
    with open(bf, "r") as f:
        content = f.read()

    # Split into hits (each hit starts with ">")
    hits = re.split(r"^>", content, flags=re.MULTILINE)[1:]
    file_hits = []

    for hit in hits:
        header = hit.splitlines()[0]
        score_match = re.search(r"Score\s+=\s+([\d.]+)\s+bits", hit)
        if score_match:
            score = float(score_match.group(1))
            file_hits.append((score, header))

    # Sort descending by score
    file_hits.sort(reverse=True, key=lambda x: x[0])

    # Take top N from this file
    result.extend(file_hits[:top_hits_per_file])

# Show all selected hits
print("Selected hits (top 2 from each file):")
for score, header in result:
    print(f"{score} bits -> {header}")
