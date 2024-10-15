import sys
from Bio import SeqIO
from Bio import SearchIO

def read_uniprot_ids(query_fasta):
    print(f"Reading Uniprot IDs from {query_fasta}...")
    uniprot_ids = set()
    for record in SeqIO.parse(query_fasta, "fasta"):
        if '|' in record.description:
            uniprot_ids.add(record.description.split("|")[1])
    print(f"Found {len(uniprot_ids)} Uniprot IDs.")
    return uniprot_ids

def filter_hits(hmmsearch_output, uniprot_ids, output_file):
    print(f"Filtering hits from {hmmsearch_output}...")
    removed_hits_ids = []
    removed_hits_count = 0
    for result in SearchIO.parse(hmmsearch_output, 'hmmer3-text'):
        for hit in result.hits:
            uniprot_id = hit.id.split('|')[1]
            if uniprot_id not in uniprot_ids:
                pass
            else:
                removed_hits_count += 1
                removed_hits_ids.append(hit.id)
    print(f"Removed {removed_hits_count} hits.")
    with open(output_file, 'w') as f:
        for hit_id in removed_hits_ids:
            f.write(hit_id + '\n')

if __name__ == "__main__":
    hmmsearch_output = sys.argv[1]
    query_fasta = sys.argv[2]
    output_file = sys.argv[3]

    print("Starting Uniprot hits removal process...")
    uniprot_ids = read_uniprot_ids(query_fasta)
    filter_hits(hmmsearch_output, uniprot_ids, output_file)
    print(f"Uniprot hits removal completed. Results written to {output_file}")
