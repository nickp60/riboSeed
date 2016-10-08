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

import logging
logger = logging

from pyutilsnrw.utils3_5 import make_output_prefix, check_installed_tools,\
    copy_file, get_ave_read_len_from_fastq, get_number_mapped,\
    extract_mapped_and_mappedmates, keep_only_first_contig, md5,\
    combine_contigs, clean_temp_dir, md5

from riboSeed.riboSeed import  check_smalt_full_install,\
    map_to_ref_smalt, convert_bams_to_fastq, estimate_distances_smalt



@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class utils3_5TestCase(unittest.TestCase):
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                        "output_riboseed_tests")
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.ref_fasta = os.path.join(self.ref_dir,
                                       'cluster1.fasta')
        self.ref_Ffastq = os.path.join(self.ref_dir,
                                        'toy_reads1.fq')
        self.ref_Rfastq = os.path.join(self.ref_dir,
                                        'toy_reads2.fq')
        self.smalt_exe = "smalt"
        self.test_estimation_file = os.path.join(self.test_dir,
                                                 "est_distance.sam")
        self.map_results_prefix = os.path.join(self.test_dir,
                                               "test_mapping")
        # self.ref_disctance
        # self.pileup
        pass

    def test_make_testing_dir(self):
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir)
        self.assertTrue(os.path.exists(self.test_dir))

    # def test_references_md5(self):
    #     """ is this paranoia, as well as bad testing?
    #     """
    #     test_pairs = [["3ba332f8a3b5d935ea6c4e410ccdf44b",
    #                    "references/combined_contigs_reference.fa"]]
    #     for i in test_pairs:
    #         self.assertEqual(i[0],
    #                          md5(os.path.join(os.path.dirname(__file__),
    #                                           i[1])))

    def test_estimate_distances_smalt(self):
        # test_estimation_file = os.path.join(self.test_dir,
        #                                     "est_distance.sam")
        if os.path.exists(self.test_estimation_file):
            print("warning! existing distance esimation file!")

        est_file = estimate_distances_smalt(outfile=self.test_estimation_file,
                                            smalt_exe=self.smalt_exe,
                                            ref_genome=self.ref_fasta,
                                            fastq1=self.ref_Ffastq,
                                            fastq2=self.ref_Rfastq,
                                            cores=1,
                                            logger=None)

    def test_map_to_ref_smalt(self):
        est_file = estimate_distances_smalt(outfile=self.test_estimation_file,
                                            smalt_exe=self.smalt_exe,
                                            ref_genome=self.ref_fasta,
                                            fastq1=self.ref_Ffastq,
                                            fastq2=self.ref_Rfastq,
                                            cores=1,
                                            logger=None)
        map_to_ref_smalt(ref=self.ref_fasta,
                         fastq_read1=self.ref_Ffastq,
                         fastq_read2=self.ref_Rfastq,
                         distance_results=est_file,
                         map_results_prefix=self.map_results_prefix,
                         cores=1, samtools_exe="samtools",
                         smalt_exe="smalt", fastq_readS="",
                         read_len=100, step=3, k=5,
                         scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
                         logger=logger)
        ####
        #TODO make actual tests;  best way to do that?


    # def test_convert_bams_to_fastq(self):
    #     convert_bams_to_fastq(map_results_prefix,
    #                           fastq_results_prefix,
    #                           keep_unmapped)

    # def test_run_spades(self):
    #     run_spades(output, ref, ref_as_contig, pe1_1='', pe1_2='', pe1_s='',
    #                as_paired=True, keep_best=True, prelim=False,
    #                groom_contigs='keep_first',
    #                k="21,33,55,77,99", seqname='', spades_exe="spades.py")

    # def test_check_samtools_pileup(self):
    #     check_samtools_pileup(self.pileup)

    # def test_reconstruct_seq(self):
    #     reconstruct_seq(refpath, pileup, verbose=True, veryverb=False,
    #                 logger=None)

    def tearDown(self):
        # test_estimation_file = os.path.join(self.test_dir,
        #                                     "est_distance.sam")
        if os.path.exists(self.test_estimation_file):
            print("removing test distance estimation file")
            os.remove(self.test_estimation_file)
        pass

if __name__ == '__main__':
    samtools_exe = "samtools"
    smalt_exe = "smalt"
    if not check_installed_tools(smalt_exe, hard=False):
        print("Error! Cannot test without install of SMALT in PATH")
    logger=logging
    unittest.main()
