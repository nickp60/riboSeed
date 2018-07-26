# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import os
import unittest
import time

# I hate this line but it works :(
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))

from riboSeed.riboSelect import get_loci_list_for_features, \
    dict_from_jenks, count_feature_hits_per_sequence, parse_args_clusters
from riboSeed.riboSnag import Locus

logger = logging
sys.dont_write_bytecode = True


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class riboSelect_TestCase(unittest.TestCase):
    """ still need to write tests for  count_feature_hits_per_sequence
    """
    def setUp(self):
        self.curdir = os.getcwd()
        self.testdirname = os.path.join(os.path.dirname(__file__),
                                        "output_riboSelect_tests")
        self.test_loci_file = os.path.join(os.path.dirname(__file__),
                                           str("references" + os.path.sep +
                                               'grouped_loci_reference.txt'))
        self.test_gb_file = os.path.join(os.path.dirname(__file__),
                                         str("references" + os.path.sep +
                                             'NC_011751.1.gb'))
        self.test_loci_file = os.path.join(os.path.dirname(__file__),
                                           str("references" + os.path.sep +
                                               'grouped_loci_reference.txt'))
        self.startTime = time.time()
        if not os.path.exists(self.testdirname):
            os.makedirs(self.testdirname, exist_ok=True)

    def test_count_feature_hits_per_sequence_all_features(self):
        nfeat_simple = count_feature_hits_per_sequence(
            all_feature=True,
            specific_features=[],
            rec_id_list=[],
            loci_list=[],
            logger=logger)
        self.assertEqual(nfeat_simple, None)

    def test_count_feature_hits_per_sequence(self):
        filtered_loci_list = [
            Locus(end_coord=228162,
                  strand= None,
                  product= '16S',
                  index= 0,
                  start_coord= 226621,
                  sequence_id='NC_011751.1',
                  locus_tag='ECUMN_16S_3',
                  feature_type='rRNA'),
            Locus(end_coord= 231513, strand= None, product= '23S', index= 1, start_coord= 228609, sequence_id='NC_011751.1', locus_tag='ECUMN_23S_3', feature_type='rRNA'),
            Locus(end_coord= 3022180, strand= None, product= '23S', index= 2, start_coord= 3019276, sequence_id='NC_011751.1', locus_tag='ECUMN_23S_2', feature_type='rRNA'),
            Locus(end_coord= 3024067, strand= None, product= '16S', index= 3, start_coord= 3022526, sequence_id='NC_011751.1', locus_tag='ECUMN_16S_2', feature_type='rRNA'),
            Locus(end_coord= 3873035, strand= None, product= '23S', index= 4, start_coord= 3870131, sequence_id='NC_011751.1', locus_tag='ECUMN_23S_1', feature_type='rRNA'),
            Locus(end_coord= 3875023, strand= None, product= '16S', index= 5, start_coord= 3873482, sequence_id='NC_011751.1', locus_tag='ECUMN_16S_1', feature_type='rRNA'),
            Locus(end_coord= 4430216, strand= None, product= '16S', index= 6, start_coord= 4428675, sequence_id='NC_011751.1', locus_tag='ECUMN_16S_4', feature_type='rRNA'),
            Locus(end_coord= 4433474, strand= None, product= '23S', index= 7, start_coord= 4430571, sequence_id='NC_011751.1', locus_tag='ECUMN_23S_4', feature_type='rRNA'),
            Locus(end_coord= 4523730, strand= None, product= '16S', index= 8, start_coord= 4522189, sequence_id='NC_011751.1', locus_tag='ECUMN_16S_5', feature_type='rRNA'),
            # Locus(end_coord= 4527078, strand= None, product= '23S', index= 9, start_coord= 4524177, sequence_id='NC_011751.1', locus_tag='ECUMN_23S_5', feature_type='rRNA'),
            Locus(end_coord= 4657586, strand= None, product= '16S', index= 10, start_coord= 4656045, sequence_id='NC_011751.1', locus_tag='ECUMN_16S_6', feature_type='rRNA'),
            Locus(end_coord= 4660845, strand= None, product= '23S', index= 11, start_coord= 4657941, sequence_id='NC_011751.1', locus_tag='ECUMN_23S_6', feature_type='rRNA'),
            Locus(end_coord= 4698680, strand= None, product= '16S', index= 12, start_coord= 4697139, sequence_id='NC_011751.1', locus_tag='ECUMN_16S_7', feature_type='rRNA'),
            Locus(end_coord= 4701939, strand= None, product= '23S', index= 13, start_coord= 4699035, sequence_id='NC_011751.1', locus_tag='ECUMN_23S_7', feature_type='rRNA')]


        nfeat_simple2 = count_feature_hits_per_sequence(
            all_feature=False,
            specific_features=["16S", "23S"],
            loci_list=filtered_loci_list,
            rec_id_list = ['NC_011751.1'],
            logger=logger)
        self.assertEqual(nfeat_simple2['NC_011751.1'], [7, 6])

    def test_dict_from_jenks(self):
        res_dict = {
            "1": [5001, 6988, 9985],
            "2": [20096, 21991, 24988],
            "3": [35099, 37086, 40083, 40328],
            "4": [50439, 52334, 55333],
            "5": [65444, 67431, 70430],
            "6": [80541, 82436, 85435],
            "7": [95546, 97441, 100440]
        }

        coli_data  = [5001, 6988, 9985, 20096, 21991, 24988, 35099, 37086,
                      40083, 40328, 50439, 52334, 55333, 65444, 67431, 70430,
                      80541, 82436, 85435, 95546, 97441, 100440]
        self.assertEqual(res_dict,
                         dict_from_jenks(
                             data=coli_data, centers=7, logger=logger))

    def test_dict_from_jenks_1cents(self):

        coli_data  = [5001, 6988, 9985, 20096, 21991, 24988, 35099, 37086,
                      40083, 40328, 50439, 52334, 55333, 65444, 67431, 70430,
                      80541, 82436, 85435, 95546, 97441, 100440]
        res_dict  = {
            "1": [5001, 6988, 9985, 20096, 21991, 24988, 35099, 37086,
                  40083, 40328, 50439, 52334, 55333, 65444, 67431, 70430,
                  80541, 82436, 85435, 95546, 97441, 100440]
            }
        self.assertEqual(res_dict,
                         dict_from_jenks(
                             data=coli_data, centers=1, logger=logger))

    def test_dict_from_jenks_0cents(self):
        coli_data  = [5001, 6988, 9985]
        with self.assertRaises(AssertionError):
            dict_from_jenks(data=coli_data, centers=0, logger=logger)

    def test_dict_from_jenks_exceed(self):

        coli_data  = [5001, 6988, 9985]
        res_dict = {
            "1": [5001],
            "2": [6988],
            "3": [9985]
        }
        self.assertEqual(
            res_dict,
            dict_from_jenks(data=coli_data, centers=5, logger=logger))

    def test_locus_tag_dict(self):
        """test different scenarios of trying to get features
        from a gb file with get_loci_list_for_features
        """
        # with self.assertRaises(SystemExit):
        # test that when given a bad feature empty results are returned
        bum_coords, bum_nfeat_simple = \
            get_loci_list_for_features(gb_path=self.test_gb_file,
                                        feature="RRNA",
                                        specific_features="16S",
                                        verbose=False,
                                        logger=logger)
        self.assertTrue([len(x) == 0 for x in
                         [bum_nfeat_simple, bum_coords]])
        bum_coords2,  bum_nfeat_simple2 = \
            get_loci_list_for_features(gb_path=self.test_gb_file,
                                        feature="rRNA",
                                        specific_features="18S",
                                        verbose=False,
                                        logger=logger)
        self.assertTrue([len(x) == 0 for x in
                         [bum_nfeat_simple2, bum_coords2]])

        filtered, nfeat_simple = \
            get_loci_list_for_features(gb_path=self.test_gb_file,
                                        feature="rRNA",
                                        specific_features="16S",
                                        verbose=True,
                                        logger=logger)
        self.assertEqual(len(filtered), 7)
        sample_entry = [x for x in filtered if \
                        x.sequence_id =='NC_011751.1' and \
                        x.start_coord == 4428675][0]
        self.assertEqual(sample_entry.locus_tag, 'ECUMN_16S_4')
        self.assertEqual(sample_entry.feature_type, 'rRNA')
        self.assertEqual(sample_entry.product, "16S")
        self.assertEqual(nfeat_simple, {'NC_011751.1': [7]})

    def test_parse_clusters_infer(self):
        centers = parse_args_clusters(clusters='', nrecs=1,
                                      logger=logger)
        self.assertEqual(centers, [0])

    def test_parse_clusters_provide(self):
        centers2 = parse_args_clusters(clusters='3', nrecs=1,
                                       logger=logger)
        self.assertEqual(centers2, [3])

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))


if __name__ == '__main__':
    unittest.main()
