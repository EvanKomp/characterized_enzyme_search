#!/bin/bash

# Check if all required arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_hmm> <output_prefix> <database>"
    exit 1
fi

# Assign input arguments to variables
inname="$1"
outname="$2"
dbname="$3"
E="$4"

# Run hmmsearch
hmmsearch \
    --cpu 40 \
    -E $E \
    --incE $E \
    --tblout "${outname}.tbl" \
    --domtblout "${outname}.domtbl" \
    --noali \
    --acc \
    --notextw \
    "${inname}" \
    "${dbname}" > "${outname}.out"

