#!/bin/bash

# Extract headers from bacteria file
grep ">" bacteria.231.1.genomic.fna > bacteria.231.1.genomic.fna.names

echo "First 10 headers:"
head bacteria.231.1.genomic.fna.names
echo "----------------------------------------------------------------------------"

echo "Last 10 headers:"
tail bacteria.231.1.genomic.fna.names
echo "-----------------------------------------------------------------------------"

echo "Total number of sequences:"
wc -l bacteria.231.1.genomic.fna.names
echo "------------------------------------------------------------------------------"
