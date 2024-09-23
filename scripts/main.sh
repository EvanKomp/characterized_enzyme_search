#!/bin/bash

# Source the parse_yaml function
source parse_yaml.sh

# Read parameters
eval $(parse_yaml params.yaml)

# Generate base name from input alignment
base_name=$(basename "$input_alignment" | sed 's/\.[^.]*$//')
echo $base_name

# Generate output file names
hmm_output="${base_name}.hmm"
hmmsearch_output="${base_name}_search"
filtered_hits1="${base_name}_id_filter1.txt"
filtered_hits2="${base_name}_id_filter2.txt"
retrieved_sequences="${base_name}_retrieved_sequences.fasta"
final_alignment="${base_name}_hits_aligned_to_hmm.a2m"

# Run pipeline steps
#####
# build hmm
hmmbuild "$hmm_output" "$input_alignment"

# search the db for hits
./hmmsearch.sh "$hmm_output" "$hmmsearch_output" "$target_db"

# remove hits in the original alignment
python remove_uniprot_hits.py "${hmmsearch_output}.out" "$input_alignment" "$filtered_hits1"

# remove hits with poor coverage
python remove_poor_coverage.py "${hmmsearch_output}.out" "$hmm_output" "$min_coverage" "$filtered_hits2"

# get the sequences for hits
python retrieve_sequences.py "${hmmsearch_output}.out" "$filtered_hits1" "$filtered_hits2" "$target_db" "$retrieved_sequences"

# align hits to the hmm
hmmalign --trim --amino --outformat a2m "$hmm_output" "$retrieved_sequences" > "$final_alignment"
