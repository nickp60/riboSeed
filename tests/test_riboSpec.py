# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import shutil
import networkx as nx
import os
import unittest
import random
from Bio import SeqIO
from argparse import Namespace
from pyutilsnrw.utils3_5 import md5
# I hate this line but it works :(
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))


from riboSeed import riboSpec as rs

sys.dont_write_bytecode = True

logger = logging


class RiboSpecTest(unittest.TestCase):
    """
    """
    def setUp(self):
        self.results_dir = os.path.join(os.path.dirname(__file__),
                                        "output_riboSpec_tests")
        self.ref_dir = os.path.join(
            os.path.dirname(__file__), "references", "")
        self.spec_ref_dir = os.path.join(
            os.path.dirname(__file__),
            "references",
            "riboSpec_references", "")
        self.spades_dir = os.path.join(
            os.path.dirname(__file__),
            "references",
            "riboSpec_references", "spades_test", "")
        self.fastg = os.path.join(self.spec_ref_dir, "mini.fastg")
        self.to_be_removed = []

    def test_graph(self):
        """ points are direction
               1
         (1) ---. (2) .  8
          .\     /     --\
        4  \  ./ 2      (3)
            (4)
          .   \------- 8---9
      16  /
        /
      (5) --------- .(6) ---- 7
             32



        """
        DG = nx.DiGraph()
        cons = [
            (1,2, 1),
            (1,4, .25),
            (2,4, 2),
            (4,1, 4),
            (4,8,1),
            (8,9, 50),
            (3,2, 8),
            (1,5, 16),
            (5,6, 32),
            (6,7, .5)]

        for i in range(len(cons)):
            DG.add_weighted_edges_from([(cons[i][0], cons[i][1], cons[i][2])])
            DG.add_weighted_edges_from([(cons[i][1], cons[i][0], cons[i][2])])
        print( "ze graph:------------------------------------------")
        print(DG.nodes.data())
        print(DG.edges.data())
        print( "------------------------------------------")
        # print(nx.dijkstra_predecessor_and_distance(DG, 4))
        # incl, border = rs.neighborhood_by_length(DG, source=4, cutoff=10)


    def test_spades_parse(self):
        t = rs.get_fastgs_from_spades_dir(self.spades_dir)
        self.assertEqual(4, len(t.keys()))

    def test_extract_node_len_cov_rc(self):
        node_name = "EDGE_313_length_1481_cov_1436.86"
        node_name_rc = "EDGE_313_length_1481_cov_1436.86'"
        self.assertEqual("313", rs.extract_node_len_cov_rc(node_name)[0])
        self.assertEqual(1481, rs.extract_node_len_cov_rc(node_name)[1])
        self.assertEqual(1436.86, rs.extract_node_len_cov_rc(node_name)[2])
        self.assertTrue(rs.extract_node_len_cov_rc(node_name_rc)[3])
        self.assertFalse(rs.extract_node_len_cov_rc(node_name)[3])

    def test_make_Node(self):
        node_name = "EDGE_313_length_1481_cov_1436.86"
        n = rs.make_Node(node_name)
        self.assertEqual(313, n.name)
        self.assertEqual(1481, n.length)
        self.assertEqual(1436.86, n.cov)
        self.assertEqual(node_name, n.raw)


    def test_make_adjacency_matrix(self):
        d = {1: [3, 4],
             2: [],
             3: [],
             4: [2, 3]}
        M = rs.make_adjacency_matrix(d)
        ref = [[0, 0, 1, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 1, 1, 0]]
        self.assertEqual(ref, M)

    def test_parse_fastg(self):
        node_list, g, DG = rs.parse_fastg(self.fastg)
        self.assertEqual([1, 2, 3,4], list(DG.nodes()))
        # no edges are added yet
        self.assertEqual(0, DG.number_of_edges())

    def test_add_temp_edges(self):
        node_list, g, DG = rs.parse_fastg(self.fastg)
        rs.add_temp_edges(node_list, DG)
        self.assertEqual(8, DG.number_of_edges())

    def tearDown(self):
        """
        """
        for filename in self.to_be_removed:
            try:
                os.unlink(filename)
            except IsADirectoryError:
                shutil.rmtree(filename)


if __name__ == '__main__':
    unittest.main()
