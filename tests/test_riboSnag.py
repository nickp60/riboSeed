# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas
The Goal of this is to have a unified place to put the useful
python 3.5 functions or templates

how I got the fastq file
# seqtk sample -s 27 ~/GitHub/FA/pseudochromosome/data/20150803_Abram1/ \
    reads/3123-1_1_trimmed.fastq .0005

bam file was from a riboseed mapping; md5: 939fbf2c282091aec0dfa278b05e94ec

mapped bam was made from bam file with the following command
 samtools view -Bh -F 4 /home/nicholas/GitHub/FB/Ecoli_comparative_genomics/
    scripts/riboSeed_pipeline/batch_coli_unpaired/map/
    mapping_20160906_region_7_riboSnag/
    test_smalt4_20160906_region_7_riboSnagS.bam >
     ~/GitHub/pyutilsnrw/tests/test_mapped.sam
md5: 27944249bf064ba54576be83053e82b0

"""
__version__ = "0.0.3"
import time
import sys
import shutil
import logging
import subprocess
import os
import unittest
import hashlib
import glob
import argparse
sys.dont_write_bytecode = True

from riboSeed.riboSnag import parse_clustered_loci_file


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
class riboSnag_TestCase(unittest.TestCase):
    # def init(self):
    #     pass
    # keep_temps = False
    # self.samtools = "test"
    # # samtools_exe = "samtoolz"
    # # def __init__(self, testname, keep_temps):
    # #     super(utils3_5TestCase, self).__init__(testname)
    # #     self.keep_temps = keep_temps
    # #     pass
    # print(self.samtools)
    def setUp(self):
        pass

    def test_parse_loci(self):
        clusters = parse_clustered_loci_file(test_loci_file, logger=logger)
        self.assertEqual(clusters[0][0], "CM000577.1")


    # # def test_this_fails(self):
    # #      self.assertEqual("pinecone", 42)

    # def test_clean_temp_dir(self):
    #     """ I tried to do something like
    #     @unittest.skipUnless(clean_temp, "temporary files were retained")
    #     but couldnt get the variabel to be passed through.
    #     """
    #     if not os.path.exists(os.path.join(testdirname, "test_subdir")):
    #         os.makedirs(os.path.join(testdirname, "test_subdir"))
    #     clean_temp_dir(testdirname)

    # def test_make_output_prefix(self):
    #     test_prefix = make_output_prefix(testdirname, "utils_3.5")
    #     self.assertEqual(test_prefix,
    #                      "".join([testdirname, os.path.sep, "utils_3.5"]))

    # def test_check_installed_tools(self):
    #     """is pwd on all mac/linux systems?
    #     #TODO replace with better passing test
    #     """
    #     check_installed_tools(["pwd"])
    #     # test fails properly
    #     with self.assertRaises(SystemExit):
    #         check_installed_tools(["thisisnotapathtoanactualexecutable"])

    # def test_md5_strings(self):
    #     """ minimal md5 examples
    #     """
    #     self.assertEqual(md5("thisstringisidenticalto", string=True),
    #                      md5("thisstringisidenticalto", string=True))
    #     self.assertNotEqual(md5("thisstringisntidenticalto", string=True),
    #                         md5("thisstringisnotidenticalto", string=True))

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
    samtools_exe = "samtools"
    unittest.main()
