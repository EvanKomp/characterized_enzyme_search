import sys
from Bio import SeqIO
from Bio.SearchIO import parse

def read_bad_ids(file_paths):
    """
    Read and combine bad IDs from multiple files.
    
    :param file_paths: List of file paths containing bad IDs
    :return: Set of bad IDs
    """
    bad_ids = set()
    for file_path in file_paths:
        with open(file_path, 'r') as f:
            bad_ids.update(line.strip() for line in f)
    return bad_ids

def retrieve_sequences(hmmer_output, bad_id_files, database_file, output_file):
    """
    Retrieve sequences from HMMER output, excluding bad IDs.
    
    :param hmmer_output: Path to HMMER output file
    :param bad_id_files: List of files containing bad IDs
    :param database_file: Path to the sequence database file
    :param output_file: Path to the output file for retrieved sequences
    """
    print(f"Processing HMMER output: {hmmer_output}")
    print(f"Reading bad IDs from: {', '.join(bad_id_files)}")
    
    bad_ids = read_bad_ids(bad_id_files)
    print(f"Found {len(bad_ids)} bad IDs to exclude.")
    
    # Parse HMMER output and filter hits
    good_hits = set()
    for query in parse(hmmer_output, 'hmmer3-text'):
        for hit in query.hits:
            if hit.id not in bad_ids:
                good_hits.add(hit.id)
    
    print(f"Found {len(good_hits)} good hits to retrieve.")
    
    # Retrieve and write sequences
    retrieved_count = 0
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(database_file, "fasta"):
            if record.id in good_hits:
                SeqIO.write(record, out_f, "fasta")
                retrieved_count += 1
    
    print(f"Retrieved and wrote {retrieved_count} sequences to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python retrieve_sequences.py <hmmer_output> <bad_ids_file1> <bad_ids_file2> <database_file> <output_file>")
        sys.exit(1)
    
    hmmer_output = sys.argv[1]
    bad_id_files = [sys.argv[2], sys.argv[3]]
    database_file = sys.argv[4]
    output_file = sys.argv[5]

    retrieve_sequences(hmmer_output, bad_id_files, database_file, output_file)