# scripts/predict_topt.py
'''
* Author: Evan Komp
* Created: 10/16/2024
* Company: National Renewable Energy Lab, Bioeneergy Science and Technology
* License: MIT

This script relies on DeepET: https://onlinelibrary.wiley.com/doi/full/10.1002/pro.4480, download and extract it, install the env.

The tools was wrapped using AIDE into thermal prediction benchmarking pipeline: https://github.com/EvanKomp/enzyme_thermal_benchmark
Run the setup there and set environmental variables DEEPET_ENV_NAME and DEEPET_PATH.
'''
import os
from tools.deepet import DeepET
from aide_predict import ProteinSequences
import pandas as pd

import dvc.api
PARAMS = dvc.api.params_show()

def main():

    ps = ProteinSequences.from_fasta(os.path.join('outputs', 'combined_hits.fasta'))
    new_ps = []
    for p in ps:
        if len(p) < PARAMS['max_seq_len']:
            new_ps.append(p)
        else:
            print(f"Skipping {p.id} due to length {len(p)}")
    ps = ProteinSequences(new_ps)

    deepet = DeepET()
    deepet.fit([])

    topts = []
    for p_batch in ps.iter_batches(1000):
        try:
            topt = deepet.predict(p_batch)
            topts.extend(list(topt.flatten()))
        except:
            max_length = max([len(p) for p in p_batch])
            print(f"Failed to predict topt for batch of {len(p_batch)} proteins with max length {max_length}")
            raise

    topt_df = pd.DataFrame({
        'ids': ps.ids,
        'topt': topts
    })

    topt_df.to_csv(os.path.join('outputs', 'topt_predictions.csv'), index=False)

if __name__ == '__main__':
    main()


