import sys
from Bio import SearchIO

def get_hmm_length(hmm_file):
    print(f"Reading HMM length from {hmm_file}...")
    with open(hmm_file, 'r') as f:
        for line in f:
            if line.startswith('LENG'):
                length = int(line.split()[1])
                print(f"HMM length: {length}")
                return length
    raise ValueError("HMM length not found in the file")

def filter_hits(input_file, hmm_length, min_coverage, output_file):
    print(f"Filtering hits from {input_file} with minimum coverage of {min_coverage}...")
    removed_hits_ids = []
    removed_hits = 0
    with open(output_file, 'w') as outfile:
        for result in SearchIO.parse(input_file, 'hmmer3-text'):
            for hit in result.hits:
                coverage = sum(hsp.hit_span for hsp in hit.hsps) / hmm_length
                if coverage >= min_coverage:
                    pass
                else:
                    removed_hits += 1
                    removed_hits_ids.append(hit.id)
    print(f"Hits removed due to poor coverage: {removed_hits}")
    with open(output_file, 'w') as f:
        for hit_id in removed_hits_ids:
            f.write(hit_id + '\n')

if __name__ == "__main__":
    input_file = sys.argv[1]
    hmm_file = sys.argv[2]
    min_coverage = float(sys.argv[3])
    output_file = sys.argv[4]

    print("Starting poor coverage hits removal process...")
    hmm_length = get_hmm_length(hmm_file)
    filter_hits(input_file, hmm_length, min_coverage, output_file)
    print(f"Poor coverage hits removal completed. Accessions written to {output_file}")
