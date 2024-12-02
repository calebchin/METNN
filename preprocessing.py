import numpy as np

"""

Computes the empirical distributions for each sample given a biom feature table
with several samples and OTUs.

Inputs:

table: a biom feature table

Returns:

empirical: an array of arrays. Each inner array is the empirical distribution 
for a given sample. The nth inner array corresponds to the sample in the nth 
column of the feature table. Likewise, the ith index in each inner array 
corresponds to the OTU in the ith row of the biom feature table.


"""

def get_empirical (table):
  data_matrix = table.matrix_data.toarray()
  empirical = []
  for sample in range(len(data_matrix[0])):
    sum_sample = 0
    for species in range(len(data_matrix)):
      sum_sample+=data_matrix[species][sample]
    sample_est = []
    for species in range(len(data_matrix)):
      sample_est.append(data_matrix[species][sample]/sum_sample)
    empirical.append(sample_est)

    return empirical
