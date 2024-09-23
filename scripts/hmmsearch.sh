#!/bin/bash

# Check if all required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_hmm> <output_prefix> <database>"
    exit 1
fi

# Assign input arguments to variables
inname="$1"
outname="$2"
dbname="$3"

# Run hmmsearch
hmmsearch \
    --cpu 40 \
    -E 1e-10 \
    --incE 1e-10 \
    -A "${outname}.sto" \
    --tblout "${outname}.tbl" \
    --domtblout "${outname}.domtbl" \
    --noali \
    --acc \
    --notextw \
    "${inname}" \
    "${dbname}" > "${outname}.out"

