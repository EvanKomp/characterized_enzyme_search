# characterized_enzyme_search

Gather homologs from a set of characterized enzymes.

## Setup

1. Install required dependencies (hmmer, mmseqs, biopython)
2. Configure `params.yaml` with input files and parameters
3. Install DVC for reproducibility

## Pipeline

1. **HMM Building**: `hmmbuild`
   - Input: Alignment file
   - Output: `characterized.hmm`

2. **HMM Search**: `hmmsearch.sh`
   - Search HMM against target databases
   - Outputs: `.tbl`, `.domtbl`, `.out` files

3. **Hit Filtering**:
   a. `remove_hits_starting_sequences.py`: Remove input sequences
   b. `remove_poor_coverage.py`: Filter based on coverage

4. **Sequence Retrieval**: `retrieve_sequences.py`
   - Retrieve filtered sequences from databases

5. **Sequence Combination**: `combine_fasta.py`
   - Combine retrieved sequences
   - Output: `combined_hits.fasta`

6. **Multiple Sequence Alignment**: `hmmalign`
   - Align retrieved sequences to original HMM
   - Output: `combined_hits_aligned.a2m`

7. **Sequence Clustering**: `mmseqs easy-cluster`
   - Cluster retrieved sequences
   - Outputs: `clustered_sequences_*.fasta/tsv`

8. **Sequence Comparison**: `compute_pid_coverage.py`
   - Compare hits to input sequences
   - Output: `hit_comparisons.tsv`

9. **Enzyme Temperature Prediction**: `predict_topt.py`
   - Predict optimal temperature for enzymes
   - Output: `topt_predictions.csv`

## Usage

Run the pipeline using DVC:

```
dvc repro
```

## Note

- Pipeline uses DVC for reproducibility
- Parameters configurable in `params.yaml`
- Enzyme temperature prediction requires a separate environment setup