# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import os
import unittest

# I hate this line but it works :(
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))

from riboSeed.riboSelect import get_filtered_locus_tag_dict, \
    dict_from_jenks, count_feature_hits, parse_args_clusters
from pyutilsnrw.utils3_5 import get_genbank_record, check_installed_tools

logger = logging
sys.dont_write_bytecode = True


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class riboSelect_TestCase(unittest.TestCase):
    """ still need to write tests for  count_feature_hits
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
        if not os.path.exists(self.testdirname):
            os.makedirs(self.testdirname, exist_ok=True)

    def test_count_feature_hits_all_features(self):
        nfeatoccur, nfeat_simple = count_feature_hits(
            all_feature=True,
            gb_path=self.test_gb_file,
            # genome_seq_records=[],
            specific_features=[],
            locus_tag_dict={},
            logger=logger)
        self.assertEqual(nfeatoccur, None)
        self.assertEqual(nfeat_simple, None)

    def test_count_feature_hits(self):
        record = get_genbank_record(self.test_gb_file)[0]
        filtered, nfeat, nfeat_simple = \
            get_filtered_locus_tag_dict(gb_path=self.test_gb_file,
                                        feature="rRNA",
                                        nrecs=1,
                                        specific_features="16S:23S",
                                        verbose=False,
                                        logger=logger)
        filtered = {
            ('NC_011751.1', 4430571): [
                8618, 'NC_011751.1', 'ECUMN_23S_4', 'rRNA',
                ['ribosomal', 'RNA', '23S']],
            ('NC_011751.1', 4524177): [
                8804, 'NC_011751.1', 'ECUMN_23S_5', 'rRNA',
                ['ribosomal', 'RNA', '23S']],
            ('NC_011751.1', 3873482): [
                7560, 'NC_011751.1', 'ECUMN_16S_1', 'rRNA',
                ['ribosomal', 'RNA', '16S']],
            ('NC_011751.1', 3022526): [
                5883, 'NC_011751.1', 'ECUMN_16S_2', 'rRNA',
                ['ribosomal', 'RNA', '16S']],
            ('NC_011751.1', 3870131): [
                7554, 'NC_011751.1', 'ECUMN_23S_1', 'rRNA',
                ['ribosomal', 'RNA', '23S']],
            ('NC_011751.1', 4522189): [
                8798, 'NC_011751.1', 'ECUMN_16S_5', 'rRNA',
                ['ribosomal', 'RNA', '16S']],
            ('NC_011751.1', 4699035): [
                9127, 'NC_011751.1', 'ECUMN_23S_7', 'rRNA',
                ['ribosomal', 'RNA', '23S']],
            ('NC_011751.1', 3019276): [
                5879, 'NC_011751.1', 'ECUMN_23S_2', 'rRNA',
                ['ribosomal', 'RNA', '23S']],
            ('NC_011751.1', 228609): [
                406, 'NC_011751.1', 'ECUMN_23S_3', 'rRNA',
                ['ribosomal', 'RNA', '23S']],
            ('NC_011751.1', 4697139): [
                9123, 'NC_011751.1', 'ECUMN_16S_7', 'rRNA',
                ['ribosomal', 'RNA', '16S']],
            ('NC_011751.1', 4656045): [
                9043, 'NC_011751.1', 'ECUMN_16S_6', 'rRNA',
                ['ribosomal', 'RNA', '16S']],
            ('NC_011751.1', 4657941): [
                9047, 'NC_011751.1', 'ECUMN_23S_6', 'rRNA',
                ['ribosomal', 'RNA', '23S']],
            ('NC_011751.1', 4428675): [
                8614, 'NC_011751.1', 'ECUMN_16S_4', 'rRNA',
                ['ribosomal', 'RNA', '16S']],
            ('NC_011751.1', 226621): [
                400, 'NC_011751.1', 'ECUMN_16S_3', 'rRNA',
                ['ribosomal', 'RNA', '16S']]
        }
        nfeatoccur2, nfeat_simple2 = count_feature_hits(
            all_feature=False,
            gb_path=self.test_gb_file,
            # genome_seq_records=[record],
            specific_features=["16S", "23S"],
            locus_tag_dict=filtered,
            logger=logger)
        self.assertEqual(nfeat_simple2['NC_011751.1'], [7, 7])
        self.assertEqual(nfeatoccur2['NC_011751.1'], [['16S', 7], ['23S', 7]])

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
        from a gb file with get_filtered_locus_tag_dict
        """
        # test single genbank record
        record = get_genbank_record(self.test_gb_file)[0]
        # with self.assertRaises(SystemExit):
        # test that when given a bad feature empty results are returned
        bum_coords, bum_nfeat, bum_nfeat_simple = \
            get_filtered_locus_tag_dict(gb_path=self.test_gb_file,
                                        feature="RRNA",
                                        specific_features="16S",
                                        nrecs=1,
                                        verbose=False,
                                        logger=logger)
        self.assertTrue([len(x) == 0 for x in
                         [bum_nfeat, bum_nfeat_simple, bum_coords]])
        bum_coords2, bum_nfeat2, bum_nfeat_simple2 = \
            get_filtered_locus_tag_dict(gb_path=self.test_gb_file,
                                        feature="rRNA",
                                        nrecs=1,
                                        specific_features="18S",
                                        verbose=False,
                                        logger=logger)
        self.assertTrue([len(x) == 0 for x in
                         [bum_nfeat2, bum_nfeat_simple2, bum_coords2]])

        filtered, nfeat, nfeat_simple = \
            get_filtered_locus_tag_dict(gb_path=self.test_gb_file,
                                        nrecs=1,
                                        feature="rRNA",
                                        specific_features="16S",
                                        verbose=False,
                                        logger=logger)
        self.assertEqual(len(filtered), 7)
        # check length of 1st item (ie, occurances) is same
        self.assertEqual(len(nfeat['NC_011751.1']),
                         len(nfeat_simple['NC_011751.1']))
        self.assertEqual([8614, 'NC_011751.1', 'ECUMN_16S_4', 'rRNA',
                          ['ribosomal', 'RNA', '16S']],
                         filtered[('NC_011751.1', 4428675)])
        self.assertEqual(nfeat_simple, {'NC_011751.1': [7]})

    def test_parse_clusters_infer(self):
        record = get_genbank_record(self.test_gb_file)[0]
        centers = parse_args_clusters(clusters='', nrecs=1,
                                      logger=logger)
        self.assertEqual(centers, [0])

    def test_parse_clusters_provide(self):
        record = get_genbank_record(self.test_gb_file)[0]
        centers2 = parse_args_clusters(clusters='3', nrecs=1,
                                       logger=logger)
        self.assertEqual(centers2, [3])

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
