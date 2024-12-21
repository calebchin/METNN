import unittest
from collections import defaultdict
from Bio import Phylo
from io import StringIO
import preprocessing as prep
import numpy as np
from biom import Table
from collections import defaultdict

class TestGetDistributions(unittest.TestCase):
    def setUp(self):
        # Sample Newick tree
        newick = "((A:0.1,B:0.2)C:0.3,(D:0.4,E:0.5)F:0.6)G;"
        self.tree = Phylo.read(StringIO(newick), "newick")

        # Sample data
        data = np.array([[10, 20, 30],
                        [5, 15, 25],
                        [7, 17, 27],
                        [3, 13, 23]])

        # Feature (OTU) IDs
        feature_ids = ['A', 'B', 'D', 'E']

        # Sample IDs
        sample_ids = ['sample1', 'sample2', 'sample3']

        # Create BIOM table
        table = Table(data, feature_ids, sample_ids)

        # Sample BIOM feature table
        self.table = table

        # Selected samples
        self.selected_samples = [0, 1, 2 ]

    def test_get_distributions(self):
        abundances, lengths, edges = prep.get_distributions(self.tree, self.table, self.selected_samples)

        print(abundances)
        print(lengths)
        print(edges)
        # # Check the structure of the returned values
        # self.assertIsInstance(abundances, list)
        # self.assertIsInstance(lengths, defaultdict)
        # self.assertIsInstance(edges, defaultdict)

        # # Check some expected values
        # self.assertEqual(lengths[0], 0.1)
        # self.assertEqual(edges[0], ('C', 'A'))
        # self.assertEqual(abundances[0]['A'], 10)


# class TestDecomposeDisjointPaths(unittest.TestCase):
#     def test_decompose_disjoint_paths(self):
#         E_w = [(1, 2), (2, 3), (3, 6), (7, 9), (9, 12)]
#         E_j = [(1, 2), (3, 6)]
        
#         result = decompose_disjoint_paths(E_w, E_j)
        
#         # Check the expected number of disjoint paths containing at least one edge in E_j
#         self.assertEqual(result, 1)

# if __name__ == '__main__':
#     unittest.main()
if __name__ == '__main__':
    unittest.main()