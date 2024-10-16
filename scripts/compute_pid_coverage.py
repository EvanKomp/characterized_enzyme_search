# scripts/compute_pid_coverage.py
'''
* Author: Evan Komp
* Created: 10/15/2024
* Company: National Renewable Energy Lab, Bioeneergy Science and Technology
* License: MIT

Compute PID and coverage of hits to original sequences using mmseqs search.

Usage:
    python compute_pid_coverage.py <input_alignment> <combined_hits> <output_tsv>
'''

import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def remove_gaps(input_alignment, output_fasta):
    """Remove gaps from input alignment and write ungapped sequences."""
    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(input_alignment, 'fasta'):
            ungapped_seq = str(record.seq).replace('-', '')
            new_record = SeqRecord(Seq(ungapped_seq), id=record.id, description='')
            SeqIO.write(new_record, out_f, 'fasta')

def run_mmseqs_search(query, target, output_prefix):
    """Run mmseqs search and parse results."""
    subprocess.run(['mmseqs', 'createdb', query, f'{output_prefix}_queryDB'])
    subprocess.run(['mmseqs', 'createdb', target, f'{output_prefix}_targetDB'])
    subprocess.run(['mmseqs', 'search', f'{output_prefix}_queryDB', f'{output_prefix}_targetDB',
                    f'{output_prefix}_resultDB', 'tmp', '-s', '8.5', '--alignment-mode', '3', '--min-seq-id', '0.3'])
    subprocess.run(['mmseqs', 'convertalis', f'{output_prefix}_queryDB', f'{output_prefix}_targetDB',
                    f'{output_prefix}_resultDB', f'{output_prefix}_results.tsv',
                    '--format-output', 'query,target,pident,fident,nident,qlen,tlen,alnlen,mismatch,qcov,tcov'])

def parse_mmseqs_results(results_file, output_tsv):
    """Parse mmseqs results and write output TSV."""
    with open(results_file) as in_f, open(output_tsv, 'w') as out_f:
        out_f.write("hit_id\tnearest_neighbor_id\t%id\tcov_query\tcov_target\n")
        for line in in_f:
            fields = line.strip().split('\t')
            hit_id, nn_id, pid, qcov, tcov = fields[1], fields[0], fields[2], fields[9], fields[10]
            out_f.write(f"{hit_id}\t{nn_id}\t{pid}\t{qcov}\t{tcov}\n")

def main(input_alignment, combined_hits, output_tsv):
    """Main function to orchestrate the PID and coverage computation."""
    ungapped_query = 'ungapped_query.fasta'
    remove_gaps(input_alignment, ungapped_query)
    
    output_prefix = 'mmseqs_output'
    run_mmseqs_search(ungapped_query, combined_hits, output_prefix)
    parse_mmseqs_results(f'{output_prefix}_results.tsv', output_tsv)

    # Clean up temporary files
    subprocess.run(['rm', '-rf', 'tmp', ungapped_query, f'{output_prefix}_queryDB*', 
                    f'{output_prefix}_targetDB*', f'{output_prefix}_resultDB*', f'{output_prefix}_results.tsv'])
    
    

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])