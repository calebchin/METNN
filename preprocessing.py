import numpy as np
from Bio import Phylo

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
    sample_est = {}
    for species in range(len(data_matrix)):
      OTU = (table.ids('observation'))[species]
      sample_est[OTU] = (data_matrix[species][sample]/sum_sample)
    empirical.append(sample_est)


    #Convert empirical into dict, access edges to multiply by L using descendants

  return empirical
  

"""
Returns a list of all the edge lengths in the tree, along with, a list of lists
of empirical distributions of edges, one for each sample. The lists of empirical 
distributions are in the same order as the list of edge lengths.

Note that as defined in the paper, an empirical distribution of an edge is the 
sum of the empirical distributions of all nodes whose paths from the root node 
pass through said edge.
"""
  
def empirical_dist_edge_len (tree, table):

  empirical_dist = get_empirical(table)
  num_samples = len(empirical_dist)
  loaded_tree = Phylo.read(tree, "newick")

  Ls = []
  dists = []

  for _ in range(num_samples):
    dists.append([])

  for node in loaded_tree.find_clades():
    if node.branch_length is not None:
      Ls.append(node.branch_length)
      descendants = [node]+list(node.find_clades())

      for i in range(num_samples):
        tot = 0
        for OTU in descendants:
          if OTU in empirical_dist[i]:
            tot+=empirical_dist[i][OTU]         
        dists[i].append(tot)

    else:
      continue

  return Ls, dists



  #Loop through edges of phylo tree object

  #for each edge, calculate empirical distribution (=sum of other empirical distributions)
  



