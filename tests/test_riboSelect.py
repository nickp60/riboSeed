# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
__version__ = "0.0.1"
import sys
import logging
import subprocess
import os
import unittest
import argparse
sys.dont_write_bytecode = True

from riboSeed.riboSelect import multisplit, get_filtered_locus_tag_dict, \
    pure_python_kmeans


from pyutilsnrw.utils3_5 import get_genbank_record, check_installed_tools,\
    multisplit


logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class riboSelect_TestCase(unittest.TestCase):
    def setUp(self):
        self.curdir = os.getcwd()
        self.testdirname = os.path.join(os.path.dirname(__file__),
                                        "output_utils3_5_tests")
        self.test_loci_file = os.path.join(os.path.dirname(__file__),
                                           str("references" + os.path.sep +
                                               'grouped_loci_reference.txt'))
        self.test_gb_file = os.path.join(os.path.dirname(__file__),
                                         str("references" + os.path.sep +
                                             'NC_011751.1.gb'))
        self.test_loci_file = os.path.join(os.path.dirname(__file__),
                                           str("references" + os.path.sep +
                                               'grouped_loci_reference.txt'))

    def test_get_filtered_locus_tag_dict(self):
        # test single genbank record
        record = get_genbank_record(self.test_gb_file)[0]
        # with self.assertRaises(SystemExit):
        # test that when given a bad feature empty results are returned
        bum_coords, bum_nfeat, bum_nfeat_simple = \
            get_filtered_locus_tag_dict([record], feature="RRNA",
                                        specific_features="16S",
                                        verbose=False,
                                        logger=logger)
        self.assertTrue([len(x) == 0 for x in \
                         [bum_nfeat, bum_nfeat_simple, bum_coords]])
        bum_coords2, bum_nfeat2, bum_nfeat_simple2 = \
            get_filtered_locus_tag_dict([record], feature="rRNA",
                                        specific_features="18S",
                                        verbose=False,
                                        logger=logger)
        self.assertTrue([len(x) == 0 for x in \
                         [bum_nfeat, bum_nfeat_simple, bum_coords]])

        filtered, nfeat, nfeat_simple =\
                get_filtered_locus_tag_dict([record], feature="rRNA",
                                            specific_features="16S",
                                            verbose=False,
                                            logger=logger)
        self.assertEqual(len(filtered), 7)
        print(filtered)
        self.assertEqual([8614, 'NC_011751.1', 'ECUMN_16S_4', 'rRNA',
                          ['ribosomal', 'RNA', '16S']],
                         filtered[('NC_011751.1', 4428675)])
        self.assertEqual(nfeat_simple, {'NC_011751.1': [7]})
        # test scaffolded
        #TODO

    def test_dopey_kmeans_function(self):
        r_is_installed = check_installed_tools("R",
                                               hard=False, logger=logger)
        if r_is_installed:
            test_for_clustering = [4, 5, 5, 6, 10, 3, 18, 34, 44, 38]
            clusters = pure_python_kmeans(test_for_clustering, group_by=None,
                                          centers=3, DEBUG=True, kind=int)
            ref_dict = {'2': [3, 4, 5, 5, 6, 10], '1': [34, 38, 44], '3': [18]}
            self.assertEqual(ref_dict, clusters)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
