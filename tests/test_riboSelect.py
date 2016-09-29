# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
__version__ = "0.0.3"
import sys
import logging
import subprocess
import os
import unittest
import argparse
sys.dont_write_bytecode = True

from riboSeed.riboSelect import multisplit, get_filtered_locus_tag_dict, \
    pure_python_kmeans


from pyutilsnrw.utils3_5 import get_genbank_record, check_installed_tools


def get_args():
    parser = argparse.ArgumentParser(
        description="test suite for pyutilsnrw repo")
    parser.add_argument("-k", "--keep_temps", dest='keep_temps',
                        action="store_true",
                        help="set if you want to inspect the output files",
                        default=False)
    args = parser.parse_args()
    return(args)

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class riboSelect_TestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_multisplit(self):
        test_string = "look_this+is+a locus_that_is+multi-delimited"
        list_of_things = multisplit(["-", "_", "+", " "], test_string)
        test_other_string = "look_this+is+a\faillocus_that_is+multi-delimited"
        list_of_other_things = multisplit(["-", "_", "+", " "],
                                          test_other_string)
        self.assertEqual(list_of_things, ["look", "this", "is", "a", "locus",
                                          "that", "is", "multi", "delimited"])
        self.assertNotEqual(list_of_other_things, ["look", "this", "is", "a",
                                                   "locus", "that", "is",
                                                   "multi", "delimited"])

    def test_get_filtered_locus_tag_dict(self):
        # test single genbank record
        record = get_genbank_record(test_gb_file)
        # with self.assertRaises(SystemExit):
        # test that when given a bad feature empty results are returned
        bum_coords, bum_nfeat, bum_nfeat_simple = \
            get_filtered_locus_tag_dict([record], feature="RRNA",
                                        specific_features="16S",
                                        verbose=False)
        self.assertTrue([len(x) == 0 for x in \
                         [bum_nfeat, bum_nfeat_simple, bum_coords]])
        bum_coords2, bum_nfeat2, bum_nfeat_simple2 = \
            get_filtered_locus_tag_dict([record], feature="rRNA",
                                        specific_features="18S",
                                        verbose=False)
        self.assertTrue([len(x) == 0 for x in \
                         [bum_nfeat, bum_nfeat_simple, bum_coords]])

        filtered, nfeat, nfeat_simple =\
                get_filtered_locus_tag_dict([record], feature="rRNA",
                                            specific_features="16S",
                                            verbose=False)
        self.assertEqual(len(filtered), 7)
        self.assertEqual([8614, 'NC_011751.1', 'ECUMN_16S_4', 'rRNA',
                          ['ribosomal', 'RNA', '16S']],
                         filtered[8614])
        self.assertEqual(nfeat_simple, {'NC_011751.1': [7]})
        # test scaffolded
        #TODO

    def test_dopey_kmeans_function(self):
        r_is_installed = check_installed_tools(["R"],
                                               hard=False, logger=logger)
        if r_is_installed:
            test_for_clustering = [4, 5, 5, 6, 10, 3, 18, 34, 44, 38]
            clusters = pure_python_kmeans(test_for_clustering, group_by=None,
                                          centers=3, DEBUG=True)
            ref_dict = {'2': [3, 4, 5, 5, 6, 10], '1': [34, 38, 44], '3': [18]}
            self.assertEqual(ref_dict, clusters)

    def tearDown(self):
        pass


if __name__ == '__main__':
    args = get_args()
    curdir = os.getcwd()
    # samtools_exe = args.samtools_exe
    testdirname = os.path.join(os.path.dirname(__file__),
                               "output_utils3_5_tests")
    test_loci_file = os.path.join(os.path.dirname(__file__),
                                   str("references" + os.path.sep +
                                       'grouped_loci_reference.txt'))
    test_gb_file = os.path.join(os.path.dirname(__file__),
                                   str("references" + os.path.sep +
                                       'NC_011751.1.gb'))
    test_loci_file = os.path.join(os.path.dirname(__file__),
                                   str("references" + os.path.sep +
                                       'grouped_loci_reference.txt'))
    # utils3_5TestCase.keep_temps = args.keep_temps
    unittest.main()
