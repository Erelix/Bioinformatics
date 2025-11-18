#!/bin/bash
while read seq_id start end; do
    if [ $start -lt $end ]; then
        ./seqkit subseq --chr "$seq_id" -r "${start}:${end}" final_representatives.fasta
    else
        ./seqkit subseq --chr "$seq_id" -r "${end}:${start}" final_representatives.fasta | ./seqkit seq -r -p
    fi
done < coordinates.txt

chmod +x extract_sequences.sh