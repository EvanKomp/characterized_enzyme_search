# characterized_enzyme_search
Gather homologs from a set of characterized enzymes

## 1. HMM Building
- Build HMM from input alignment
- Output: `characterized.hmm`

## 2. HMM Search
- Search HMM against target databases
- Outputs: `.tbl`, `.domtbl`, `.out` files

## 3. Hit Filtering
### a. Remove UniProt Hits
- Filter out sequences from input alignment
- Output: `removed_uniprot_hits_{db}.txt`

### b. Remove Poor Coverage Hits
- Filter hits based on minimum coverage
- Output: `removed_poor_coverage_{db}.txt`

## 4. Sequence Retrieval
- Retrieve filtered sequences from databases
- Output: `retrieved_sequences_{db}.fasta`

## 5. Sequence Combination
- Combine retrieved sequences from all databases
- Output: `combined_hits.fasta`

## 6. Multiple Sequence Alignment
- Align retrieved sequences to original HMM
- Include original alignment
- Output: `final_alignment.a2m`

## Note
- Pipeline uses DVC for reproducibility
- Params configurable in `params.yaml`