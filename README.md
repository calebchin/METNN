# METNN 

Welcome!

## Installation
Please see the `requirements.txt` to install the necessary packages. Also, install QIIME for preprocessing. 

## Experiments instructions
There are two ways that we ran our experiments for this project:
1. Through the `experiments.py` script. To replicate our experiments, download the data files under `processed_data`. These should include our processed data: `feature-table.biom` and `tree.nwk`. If they aren't already, place these files into a directory named "processed_data" within this repos directory. Then, run either the function `unifrac_experiment` (data for k and N increase analysis) or `exper` (UniFrac vs MET comparison) within the main function of `experiments.py`. You can run the script with the command `python experiments.py`.
2. For visualizations and the kNN work, please visit `METNN_experiments.ipynb`. Apologies for the clutter! It was recently in construction. Under "NN Tests - Vary N, k", you can find the cells used to run the experiments with kNN with increasing N and k. Under `NN Tests`, you can find the cells used to test UniFrac andn MET under kNN. The cells below this were used to generate the visualizations that you see in the paper. 

## Preprocessing instructions
1. Download the data from the [IBDMDB database](https://ibdmdb.org/results). Specifically, click on the rawfiles link for HMP2, data type 16S at the top of the table. You will have to individually click on all of these links to download the data we used.
2. Next, place the files into a folder and run the `extract_manifest.py` script. This will create a manifest for the data. Next, download the SEPP 16s rRNA Greengenes 13_8 reference from QIIME at the bottom of the page[(https://docs.qiime2.org/2021.11/data-resources/)].
3. Run the following commands using QIIME

`qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path hmp2_16s_manifest.tsv \
--output-path se_hmp2_16s.qza \
--input-format SingleEndFastqManifestPhred33V2
`

`qiime dada2 denoise-single \
  --i-demultiplexed-seqs se_hmp2_16s.qza \     
  --p-trim-left 0 \ 
  --p-trunc-len 250 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
  `
`
qiime fragment-insertion sepp \
--i-representative-sequences rep-seqs.qza \
--i-reference-database sepp-refs-gg-13-8.qza \
--o-tree insertion_tree.qza \
--o-placements tree_placements.qza \
--p-threads 12
`

`qiime tools export --input-path insertion_tree.qza --output-path final_exported_tree 
`

`qiime tools export --input-path table.qza --output-path exported-feature-table 
`
