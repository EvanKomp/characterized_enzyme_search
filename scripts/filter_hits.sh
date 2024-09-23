#!/bin/bash

# Check if all required arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_sto> <query_fasta>"
    exit 1
fi

input_sto="$1"
query_fasta="$2"
output_base="${input_sto%.*}"

# Convert STO to A2M
echo "Converting STO to A2M..."
esl-reformat a2m "$input_sto" > "${output_base}.a2m"

if [ $? -ne 0 ]; then
    echo "Error: Failed to convert STO to A2M."
    exit 1
fi

# Extract query IDs
echo "Extracting query IDs..."
query_ids=$(grep "^>" "$query_fasta" | sed 's/^>//; s/|.*$//' | sort -u)

# Filter A2M file
echo "Filtering A2M file..."
temp_file="${output_base}_temp.a2m"
current_seq=""
include_seq=true

while IFS= read -r line; do
    if [[ "$line" =~ ^\> ]]; then
        if [ "$include_seq" = true ] && [ -n "$current_seq" ]; then
            echo "$current_seq"
        fi
        current_seq="$line"
        seq_id=$(echo "$line" | sed 's/^>//; s/|.*$//')
        if echo "$query_ids" | grep -q "^$seq_id$"; then
            include_seq=false
        else
            include_seq=true
        fi
    elif [ "$include_seq" = true ]; then
        current_seq+=$'\n'"$line"
    fi
done < "${output_base}.a2m" > "$temp_file"

# Output the last sequence if it wasn't a query
if [ "$include_seq" = true ] && [ -n "$current_seq" ]; then
    echo "$current_seq" >> "$temp_file"
fi

mv "$temp_file" "${output_base}_filtered.a2m"

echo "Filtered A2M file saved as ${output_base}_filtered.a2m"

# Count sequences
total_seqs=$(grep -c "^>" "${output_base}.a2m")
filtered_seqs=$(grep -c "^>" "${output_base}_filtered.a2m")
removed_seqs=$((total_seqs - filtered_seqs))

echo "Summary:"
echo "  Total sequences in original A2M: $total_seqs"
echo "  Sequences after filtering: $filtered_seqs"
echo "  Sequences removed (queries): $removed_seqs"
