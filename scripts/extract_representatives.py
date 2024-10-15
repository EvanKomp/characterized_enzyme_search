# scripts/extract_representatives.py
'''
* Author: Evan Komp
* Created: 10/9/2024
* Company: National Renewable Energy Lab, Bioeneergy Science and Technology
* License: MIT
'''
import sys
import json
import csv
import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2

def seqid(seq1, seq2):
    """Compute sequence identity between two sequences."""
    alignment = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)[0]
    return sum(a == b for a, b in zip(alignment[0], alignment[1])) / len(alignment[0])

def extract_uniprot_id(header):
    """Extract UniProt ID from FASTA header."""
    return header.split('|')[1] if '|' in header else header.split()[0]

def process_clusters(cluster_file, combined_hits, original_alignment, output_csv, output_json):
    """Process cluster file, extract representatives, and compute metrics."""
    clusters = {}
    with open(cluster_file, 'r') as f:
        for line in f:
            rep, member = line.strip().split('\t')
            if rep not in clusters:
                clusters[rep] = []
            clusters[rep].append(member)

    hit_seqs = {extract_uniprot_id(record.id): record for record in SeqIO.parse(combined_hits, 'fasta')}
    original_seqs = {record.description: record for record in SeqIO.parse(original_alignment, 'fasta')}
    
    representatives = []
    total_hits = 0
    max_cluster_size = 0
    
    with open(output_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['id', 'sequence', 'nearest_original', 'percent_identity'])
        
        for rep, members in tqdm.tqdm(clusters.items()):
            cluster_size = len(members) + 1  # Include representative
            total_hits += cluster_size
            max_cluster_size = max(max_cluster_size, cluster_size)
            rep_seq = hit_seqs[rep]
            max_seqid = 0
            nearest_orig = None
            for orig_id, orig_seq in original_seqs.items():
                curr_seqid = seqid(rep_seq.seq, orig_seq.seq)
                if curr_seqid > max_seqid:
                    max_seqid = curr_seqid
                    nearest_orig = orig_id
            
            csvwriter.writerow([rep_seq.description, str(rep_seq.seq), nearest_orig, f"{max_seqid:.2f}"])
            representatives.append(rep_seq)

    metrics = {
        "total_selected_sequences": len(representatives),
        "total_hits": total_hits,
        "max_cluster_size": max_cluster_size,
        "average_cluster_size": total_hits / len(representatives)
    }
    
    with open(output_json, 'w') as f:
        json.dump(metrics, f, indent=2)

if __name__ == "__main__":
    cluster_file = sys.argv[1]
    combined_hits = sys.argv[2]
    original_alignment = sys.argv[3]
    output_csv = sys.argv[4]
    output_json = sys.argv[5]
    process_clusters(cluster_file, combined_hits, original_alignment, output_csv, output_json)