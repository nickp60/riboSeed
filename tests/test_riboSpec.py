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
import networkx
from Bio import SeqIO
from argparse import Namespace
# I hate this line but it works :(
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))

from riboSeed.shared_methods import md5


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
        self.gff = os.path.join(self.spec_ref_dir, "NC_011751.1.gff")
        self.spec_gff = os.path.join(self.spec_ref_dir, "partial_barrnapped.gff")
        # from GAGE aureus
        self.real_fastg = os.path.join(self.spec_ref_dir, "assembly_graph.fastg")
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
        # print( "ze graph:------------------------------------------")
        # print(DG.nodes.data())
        # print(DG.edges.data())
        # print( "------------------------------------------")
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


    def test_neighborhood_by_length(self):
        node_list, g, DG = rs.parse_fastg(self.fastg)
        rs.add_temp_edges(node_list, DG)
        interior, border = rs.neighborhood_by_length(G=DG, source=2, cutoff=15, ignored_nodes=[])
        self.assertEqual(interior, {2, 1})
        self.assertEqual(border, {3, 4})

    def test_make_gff_list(self):
        gff_list = rs.make_gff_list(self.gff)
        self.assertEqual(22, len(gff_list))
        ref_line_1 = [
            "gi|218703261|ref|NC_011751.1|",
            "barrnap:0.7",
            "rRNA",
	    "226624",
	    "228161",
	    "0",
	    "+",
	    ".",
	    "Name=16S_rRNA;product=16S ribosomal RNA"]

        self.assertEqual(ref_line_1, gff_list[0])


    def test_find_rRNA_from_gffs(self):
        gff_list = rs.make_gff_list(self.spec_gff)
        a16, a23, a5 = rs.find_rRNA_from_gffs(gff_list, partial=False, logger=logger)
        self.assertEqual([1162], a16)
        self.assertEqual([280], a23)
        self.assertEqual([650], a5)

    def test_get_depth_of_big_nodes(self):
        node_list, g, DG = rs.parse_fastg(self.fastg)
        depths, ave, quarts = rs.get_depth_of_big_nodes(DG, threshold=0)
        ref_depths = [4.0, 3.0, 2.0, 1.0]
        ref_weighted_mean = 2.0
        ref_quarts = [1.75, 2.5, 3.25]
        self.assertEqual(ref_depths, depths)
        self.assertEqual(ref_weighted_mean, ave)
        self.assertEqual(ref_quarts, quarts)
        # nd = dict(DG.nodes(data=True))
        # ds = []
        # ls = []
        # for k, v in nd.items():
        #     ds.append(v['cov'])
        #     ls.append(v['length'])
        # prods = []
        # for i, l in enumerate(
        # weighted_mean = sum()/sum(ls)

    def test_populate_subgraph_from_source(self):
        # node 5 is the 16S

        node_list, g, DG = rs.parse_fastg(self.real_fastg)
        rs.add_temp_edges(node_list, DG)
        dg = networkx.DiGraph()
        init_node = node_list[8]
        interior, border = rs.neighborhood_by_length(G=DG, source=init_node.name,
                                                     cutoff=1000, ignored_nodes=[])
        # print("\n".join([str(x) for x in node_list]))
        # print(interior)
        # print(border)
        # print(init_node)
        dg.add_node(init_node.name, cov=init_node.cov,
                    length=init_node.length,
                    reverse_complimented= init_node.reverse_complimented,
                    raw=init_node.raw)
        rs.populate_subgraph_from_source(
            g=dg, root=init_node, length=0,
            cutoff=1000, node_list=node_list, counter=0, debug=True)
        # we cant direcly diff graphs
        print(dg.nodes())
        ref_nodes  = [5, 7, 15, 16, 17, 18, 19, 20, 30, 31, 40]
        self.assertEqual(sorted([x for x in dg.nodes()]), sorted(ref_nodes))
        # self.assertEqual(sorted([x for x in dg.nodes()]), sorted(DG.nodes()))

    def test_reverse_populate_subgraph_from_source(self):
        node_list, g, DG = rs.parse_fastg(self.fastg)
        dg = networkx.DiGraph()
        init_node = node_list[0]
        dg.add_node(init_node.name, cov=init_node.cov,
                    length=init_node.length, raw=init_node.raw)
        rs.reverse_populate_subgraph_from_source(
            g=dg, root=node_list[0], node_list=node_list, counter=0, debug=True)
        # THis is different from the forward version, because of the directionality of the graph
        self.assertEqual(sorted([x for x in dg.nodes()]), sorted([x for x in DG.nodes() if x in [4,1]]))

    def test_make_rRNAs_dict(self):
        gff_list = rs.make_gff_list(self.spec_gff)
        gff_list_partial = rs.make_gff_list(self.spec_gff)
        rrnas = rs.make_rRNAs_dict(gff_list, gff_list_partial)
        ref_rrnas = {
            '23S': {
                'partial': [911, 902, 431],
                'solid': [280]
            }, '5S': {
                'partial': [280, 281, 282],
                'solid': [650]
            }, '16S': {
                'partial': [431],
                'solid': [1162]
            }
        }
        self.assertEqual(ref_rrnas, rrnas)

    def test_check_rrnas_dict(self):
        gff_list = rs.make_gff_list(self.spec_gff)
        gff_list_partial = rs.make_gff_list(self.spec_gff)
        rrnas = rs.make_rRNAs_dict(gff_list, gff_list_partial)
        self.assertTrue(all(rs.check_rrnas_dict(rrnas, logger=logger)))

    def test_check_rrnas_dict_no16s(self):
        rrnas = {
            "16S": {
                "partial": [1, 2],
                "solid": []
            },
            "23S": {
                "partial": [1,2],
                "solid": []
            },
            "5S": {
                "partial": [1,2],
                "solid": [1,2]
            }
        }
        run16s, run23s = rs.check_rrnas_dict(rrnas, logger=logger)
        self.assertFalse(run16s)
        self.assertFalse(run23s)

    def test_remove_nested_lists(self):
        test_list = [
            [1,2,3,4,5],
            [1,2,3],
            [1,2,3,4,5],
            [1,5,6,7,8]
        ]
        dedup = rs.remove_duplicate_nested_lists(test_list)
        self.assertEqual(3, len(dedup))

    def test_find_collapsable_partial_rRNA_nodes(self):
        node_list, g, DG = rs.parse_fastg(self.real_fastg)
        nodes_data = dict(DG.nodes(data=True))
        rrnas = {'16S': {'partial': [], 'solid': [5]}, '23S': {'partial': [], 'solid': [1]}, '5S': {'partial': [1, 11, 37, 38, 39], 'solid': []}}
        # {
        #     '16S': {'partial': [682], 'solid': [63]},
        #     '5S': {'partial': [356, 71, 142, 914, 731, 697, 282, 283,
        #                        605, 606, 735],
        #            'solid': []},
        #     '23S': {'partial': [914, 44, 717, 718], 'solid': []}
        # }
        # these were verified visually with bandage
        collapsed = rs.find_collapsable_partial_rRNA_nodes(
            rrnas, DG,
            nodes_data=nodes_data, logger=logger)
        print(collapsed)

    def test_remove_similar_lists(self):
        """ verify list simplification works
        removing similar lists is when we look at the paths containing medium length nodes ( in this case, 200bp.  list "l" is the lengths of the nodes, and list "n" is the list of node names.  We are hoping to retain all the paths to node "5" that dont have duplicate lengths.  The magic is that we simplify the lights by removing values > threshold.
        """
        l = [
            [500,300,100,100,500],
            [500,300,100,150,500], # this should be removed
            [500,300,400,900,500],
            [100,100,400,900,800],
        ]
        n = [
            [1,2,3,4,5],
            [1,3,4,6,5],
            [7,8,6,1,5],
            [17,28,6,1,35],
        ]
        # min threshold is 200
        ref_paths = [
            [17, 28, 6, 1, 35],
            [7, 8, 6, 1, 5],
            [1, 2, 3, 4, 5]]
        self.assertEqual(
            ref_paths,
            rs.remove_similar_lists(n, l)
        )


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
