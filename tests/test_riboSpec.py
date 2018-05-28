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
