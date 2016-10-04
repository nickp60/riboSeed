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

from pyutilsnrw.utils3_5 import get_genbank_seq, get_genbank_record

from riboSeed.riboSnag import parse_clustered_loci_file, \
    extract_coords_from_locus, get_genbank_seq_matching_id,\
    stitch_together_target_regions


logger = logging
@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class riboSnag_TestCase(unittest.TestCase):
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
        self.samtools_exe = "samtools"

    def test_parse_loci(self):
        """
        """
        clusters = parse_clustered_loci_file(self.test_loci_file,
                                             logger=logger)
        self.assertEqual(clusters[0][0], "CM000577.1")

    def test_extract_coords_from_locus(self):

        records = get_genbank_record(self.test_gb_file)
        coord_list = extract_coords_from_locus(genome_seq_records=records,
                                               locus_tag_list=["ECUMN_0004"],
                                               feature="CDS",
                                               verbose=True, logger=logger)[0]
        loc_index = coord_list[0]
        locus_tag = coord_list[4]
        strand = coord_list[2]
        coords = coord_list[1]
        seqid  = coord_list[5]
        self.assertEqual(loc_index, 0)
        self.assertEqual(locus_tag, "ECUMN_0004")
        self.assertEqual(strand, 1)
        self.assertEqual(coords,[3692, 4978])
        self.assertEqual(seqid,'NC_011751.1')

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
