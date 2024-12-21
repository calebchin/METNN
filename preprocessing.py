import numpy as np
from Bio import Phylo
from collections import defaultdict

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
Computes the empirical distributions for a subset of samples given a biom feature table.

Inputs:
- table: a biom feature table
- selected_samples: list of sample indices to process

Returns:
- empirical: a list of empirical distributions for the selected samples.
"""
def get_empirical_subset(table, selected_samples):

    data_matrix = table.matrix_data.toarray()
    empirical = []

    for sample_index in selected_samples:
        sum_sample = np.sum(data_matrix[:, sample_index])
        sample_est = {}

        for species_index in range(len(data_matrix)):
            OTU = table.ids('observation')[species_index]
            sample_est[OTU] = (data_matrix[species_index][sample_index] / sum_sample)

        empirical.append(sample_est)

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

"""
Assigns default names to nodes in a Newick tree that lack explicit names.

Parameters:
- tree (Bio.Phylo.BaseTree.Tree): A parsed Newick tree.

Returns:
- Bio.Phylo.BaseTree.Tree: The tree with all unnamed nodes assigned default names.
"""
def assign_default_names(tree):
 
  unnamed_counter = 1  # Counter for generating default names

  # Traverse all nodes (clades) in the tree
  for clade in tree.find_clades():
      if clade.name is None:
          clade.name = f"Inner{unnamed_counter}"  # Assign a default name
          unnamed_counter += 1

  return tree


"""
Returns P_e for all edges in a given tree, the branch length of each edge, and a mapping from edge index to the paired node representation

Assumes input tree is a Bio.Phylo.BaseTree.Tree object with all named nodes

Parameters:
- tree (Bio.Phylo.BaseTree.Tree): A parsed Newick tree.
- table: BIOM feature table
- selected_samples: List of sample indices into the feature table to process
"""
def get_distributions(tree, table, selected_samples):
  empirical_dist = get_empirical_subset(table, selected_samples)
  num_samples = len(empirical_dist)
  #print(num_samples)
  # abundance of an edge is the sum of the abundances of all nodes that pass through it
  abundances = []
  for i in range(num_samples):
    abundances.append(defaultdict(float))
  #abundances = [defaultdict(lambda: 0.0)] * num_samples # dict indexed by edge index
  lengths = defaultdict(float) # indexed by edge index
  edges = defaultdict(tuple) # indexed by edge index
  edge_counter = 0
  def abund(clade, parent_name):
    nonlocal edge_counter
    # note: internal nodes do not have abundances
    curr_edge = edge_counter
    edge_counter += 1
    if clade.branch_length is not None:
      lengths[curr_edge] = clade.branch_length
      edges[curr_edge] = (parent_name, clade.name)
    if clade.is_terminal():
      #lengths[num_edge] = clade.branch_length
      for i in range(num_samples):
        if clade.name in empirical_dist[i]:
          abundances[i][curr_edge] = empirical_dist[i].get(clade.name, 0.0)
      return
    else:
      #print(clade.clades)
      for child in clade.clades:
        child_edge = edge_counter
        abund(child, clade.name)
        for i in range(num_samples):
          abundances[i][curr_edge] += abundances[i][child_edge]
  abund(tree.root, None)
  # in this case, 0th edge is the root edge, which has no length
  for i in range(num_samples):
    if 0 in abundances[i]:
          del abundances[i][0]
    abundances[i].default_factory = None # immutable dict
  if 0 in lengths:
      del lengths[0]
  if 0 in edges:
      del edges[0]
  lengths.default_factory = None
  edges.default_factory = None
  
  return abundances, lengths, edges



"""
Returns a subset of the empirical distributions for selected edges (OTUs) and their lengths,
based on a subset of selected samples.

Inputs:
- tree: Phylogenetic tree (in Newick format)
- table: BIOM feature table
- selected_samples: List of sample indices to process
- selected_OTUs: List of OTU IDs (edges) to process

Returns:
- Ls: List of edge lengths for the selected OTUs
- dists: List of empirical distributions for the selected edges (OTUs)
"""
def empirical_dist_edge_len_subset(tree, table, selected_samples, selected_OTUs):
  empirical_dist = get_empirical_subset(table, selected_samples)
  num_samples = len(empirical_dist)


  loaded_tree = Phylo.read(tree, "newick")

  # Initialize lists to store edge lengths and distributions
  Ls = []
  dists = [[] for _ in range(num_samples)]
  print("h")
  for node in loaded_tree.find_clades():
      if node.name in selected_OTUs:
          if node.branch_length is not None:
              Ls.append(node.branch_length)
              descendants = [node] + list(node.find_clades())
              print(node.name)
              print(list(node.find_clades()))
              for i in range(num_samples):
                  total = 0
                  for OTU in descendants:
                      if OTU.name in empirical_dist[i]:
                          total += empirical_dist[i][OTU.name]
                  dists[i].append(total)

  return Ls, dists




#Loop through edges of phylo tree object

#for each edge, calculate empirical distribution (=sum of other empirical distributions)




