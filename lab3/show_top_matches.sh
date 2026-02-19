#!/bin/bash

echo "=== Top Bacteria Matches ==="
if [ -f top_bacteria_matches.txt ]; then
    cat top_bacteria_matches.txt
else
    echo "top_bacteria_matches.txt not found"
fi

echo ""
echo "=== Top Archaea Matches ==="
if [ -f top_archaea_matches.txt ]; then
    cat top_archaea_matches.txt
else
    echo "top_archaea_matches.txt not found"
fi

echo ""
echo "=== Top Fungi Matches ==="
if [ -f top_fungi_matches.txt ]; then
    cat top_fungi_matches.txt
else
    echo "top_fungi_matches.txt not found"
fi
