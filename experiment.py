import pandas as pd
import biom
import preprocessing as prep
import MET
import UniFrac as UF
import numpy as np
from Bio import Phylo

def unifrac_experiment():
    table = biom.load_table('processed_data/feature-table.biom')
    ids = table.ids('observation')
    sampleid = table.ids('sample')
    tree = Phylo.read('processed_data/tree.nwk', "newick")
    tree = prep.assign_default_names(tree)

    df_meta = pd.read_csv('hmp2_metadata_2018-08-20.csv')
    df_meta = df_meta[df_meta['External ID'].isin(sampleid)]

    N = [10, 25, 50, 100]
    all_dists = []
    truth_labels = []
    for i, n in enumerate(N):
        random_sample = df_meta.sample(n=n, random_state=i)
        sel_samples = list(random_sample["External ID"])
        all_smpls = [np.where(sampleid == sample_id)[0][0] for sample_id in sel_samples]
        abundances, lengths, edges = prep.get_distributions(tree, table, all_smpls)
        UniFrac_dists = np.full([n, n], np.nan)
        
        counter = 0
        for i in range(len(abundances)):
            for j in range(i, len(abundances)):
                counter += 1
                print("Calculation " + str(counter) + " of 100")
            
                UniFrac_dists[i][j] = UF.weighted_unifrac(abundances[i], abundances[j], lengths)
                UniFrac_dists[j][i] = UniFrac_dists[i][j]
        truth_labels.append(list(random_sample["diagnosis"]))
        all_dists.append(UniFrac_dists)
    #np.save('MET_dists.npy', MET_dists)
    np.save('UniFrac_dists_N.npy', np.array(all_dists, dtype=object), allow_pickle=True)
    np.save("truth_labels.npy", np.array(truth_labels, dtype=object), allow_pickle=True)

def exper():

    table = biom.load_table('processed_data/feature-table.biom')
    ids = table.ids('observation')
    sampleid = table.ids('sample')
    tree = Phylo.read('processed_data/tree.nwk', "newick")
    tree = prep.assign_default_names(tree)

    df_meta = pd.read_csv('hmp2_metadata_2018-08-20.csv')
    df_meta = df_meta[df_meta['External ID'].isin(sampleid)]
    df_noibd_meta = df_meta[df_meta['diagnosis'] == 'nonIBD']
    df_ibd_meta = df_meta[df_meta['diagnosis'] == 'UC']
    
    random_sample_noibd = df_noibd_meta.sample(n=10, random_state=42)
    sel_samples_noibd = list(random_sample_noibd["External ID"])

    random_sample_ibd = df_ibd_meta.sample(n=10, random_state=42)
    sel_samples_ibd = list(random_sample_ibd["External ID"])

    # Find the corresponding indices of the selected sample IDs
    selected_sample_indices_noibd = [np.where(sampleid == sample_id)[0][0] for sample_id in sel_samples_noibd]
    selected_sample_indices_ibd = [np.where(sampleid == sample_id)[0][0] for sample_id in sel_samples_ibd]
    all_smpls = selected_sample_indices_noibd + selected_sample_indices_ibd
    # no ibd is the first 10
    abundances, lengths, edges = prep.get_distributions(tree, table, all_smpls)
    MET_dists = np.full([20, 20], np.nan)
    UniFrac_dists = np.full([20, 20], np.nan)
    counter = 0
    for i in range(len(abundances)):
        for j in range(i, len(abundances)):
            counter += 1
            print("Calculation " + str(counter) + " of 100")
            MET_dists[i][j] = MET.moment_screening_estimator(tree, abundances[i], abundances[j], lengths, 31184, edges)
            MET_dists[j][i] = MET_dists[i][j]

            UniFrac_dists[i][j] = UF.weighted_unifrac(abundances[i], abundances[j], lengths)
            UniFrac_dists[j][i] = UniFrac_dists[i][j]
    np.save('MET_dists.npy', MET_dists)
    np.save('UniFrac_dists.npy', UniFrac_dists)

def main():
    unifrac_experiment()

if __name__ == "__main__":
    main()