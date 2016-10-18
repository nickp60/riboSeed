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
    combine_contigs, clean_temp_dir, file_len

from riboSeed.riboSeed import  check_smalt_full_install,\
    map_to_ref_smalt, convert_bams_to_fastq, estimate_distances_smalt,\
    run_spades



@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class utils3_5TestCase(unittest.TestCase):
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_riboseed_tests")
        self.spades_dir = os.path.join(os.path.dirname(__file__),
                                       "output_riboseed_tests",
                                       "SPAdes_results")
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.ref_fasta = os.path.join(self.ref_dir,
                                      'cluster1.fasta')
        self.ref_Ffastq = os.path.join(self.ref_dir,
                                       'toy_reads1.fq')
        self.ref_Rfastq = os.path.join(self.ref_dir,
                                       'toy_reads2.fq')
        self.ref_bam_prefix = os.path.join(self.ref_dir,
                                       'test_bam_to_fastq')
        self.smalt_exe = "smalt"
        self.samtools_exe = "samtools"
        self.test_estimation_file = os.path.join(self.test_dir,
                                                 "est_distance.sam")
        self.map_results_prefix = os.path.join(self.test_dir,
                                               "test_mapping")
        self.fastq_results_prefix = os.path.join(self.test_dir,
                                                 "test_bam_to_fastq")
        # self.ref_disctance
        # self.pileup
        pass

    def test_make_testing_dir(self):
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir)
        self.assertTrue(os.path.exists(self.test_dir))

    def test_estimate_distances_smalt(self):
        if os.path.exists(self.test_estimation_file):
            print("warning! existing distance esimation file!")
        est_file = estimate_distances_smalt(outfile=self.test_estimation_file,
                                            smalt_exe=self.smalt_exe,
                                            ref_genome=self.ref_fasta,
                                            fastq1=self.ref_Ffastq,
                                            fastq2=self.ref_Rfastq,
                                            cores=1,
                                            logger=logger)
        ref_smi_md5 = "a444ccbcb486a8af29736028640e87cf"  # determined manually
        ref_sma_md5 = "4ce0c8b453f2bdabd73eaf8b5ee4f376"  # determined manually
        ref_mapping_len = 9271  # mapping doesnt have exact order, so cant md5
        self.assertEqual(ref_smi_md5, md5(str(est_file + ".smi")))
        self.assertEqual(ref_sma_md5, md5(str(est_file + ".sma")))
        self.assertEqual(ref_mapping_len, file_len(est_file))

    def test_convert_bams_to_fastq(self):
        if not os.path.exists(self.fastq_results_prefix):
            os.makedirs(self.fastq_results_prefix)
        convert_bams_to_fastq(map_results_prefix=self.ref_bam_prefix,
                              fastq_results_prefix=self.fastq_results_prefix,
                              keep_unmapped=False,
                              samtools_exe=self.samtools_exe)

    def test_map_to_ref_smalt(self):
        """ this tests that the resulting bam files when converted to
        sam has the right number of lines.  Not that helpful, open to
        other ideas
        """
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
                         read_len=100, step=10, k=20,
                         score_minimum=40,
                         scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
                         logger=logger)
        wc_cmds = ["samtools view -S {0} |wc".format(self.map_results_prefix +
                                                     "_pe.bam")]
        wcres = subprocess.run(wc_cmds, shell=sys.platform != "win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, check=True)
        # lines and words determinied manually.
        # Bytes change depending on cline call
        # wc returns "lines\swords\s\sbyte
        lines_and_words = ["36350", "", "472550"]
        for i in range(0, 3):
            parse_results = wcres.stdout.decode("utf-8").strip().split(" ")
            self.assertEqual(parse_results[i],
                             lines_and_words[i])

    def test_run_spades(self):
        """The tests a 'prelim' and non 'prelim' assembly against manually
        determined md5sums of the resulting contig files
        """
        contigs_ref1 ="68829130b1405e9108a02f4cd414f057"
        """PARAMS FOR CONTIGS1:
        spades.py -k 21,33,55 --trusted-contigs tests/references/cluster1.fasta
        --pe1-1 ./toy_reads1.fq --pe1-2 ./toy_reads2.fq -o spadesman --careful
        """
        contigs_ref2 = "740310315e5547e25c7aca012d198c65"
        """PARAMS FOR CONTIGS1:
        spades.py -k 21,33,55 --trusted-contigs tests/references/cluster1.fasta
        --pe1-1 ./toy_reads1.fq --pe1-2 ./toy_reads2.fq -o spadesman
        --only-assembler --cov-cutoff off --sc --careful
        manually select just the first contig NODE1..
        """
        contigs1, success = run_spades(output=os.path.join(self.spades_dir,
                                                           "test1"),
                                       ref=self.ref_fasta,
                                       ref_as_contig="trusted",
                                       pe1_1=self.ref_Ffastq,
                                       pe1_2=self.ref_Rfastq, pe1_s='',
                                       as_paired=True, keep_best=True,
                                       prelim=False,
                                       groom_contigs='keep_first',
                                       k="21,33,55", seqname='',
                                       spades_exe="spades.py", logger=logger)
        contigs2, success = run_spades(output=os.path.join(self.spades_dir,
                                                           "test2"),
                                       ref=self.ref_fasta,
                                       ref_as_contig="trusted",
                                       pe1_1=self.ref_Ffastq,
                                       pe1_2=self.ref_Rfastq, pe1_s='',
                                       as_paired=True, keep_best=True,
                                       prelim=True,
                                       groom_contigs='keep_first',
                                       k="21,33,55", seqname='',
                                       spades_exe="spades.py", logger=logger)
        self.assertEqual(contigs_ref1, md5(contigs1))
        self.assertEqual(contigs_ref2, md5(contigs2))

    # def test_check_samtools_pileup(self):
    #     check_samtools_pileup(self.pileup)

    # def test_reconstruct_seq(self):
    #     reconstruct_seq(refpath, pileup, verbose=True, veryverb=False,
    #                     logger=None)

    def tearDown(self):
        # # test_estimation_file = os.path.join(self.test_dir,
        # #                                     "est_distance.sam")
        # if os.path.exists(self.test_estimation_file):
        #     print("removing test distance estimation file")
        #     os.remove(self.test_estimation_file)
        pass

if __name__ == '__main__':
    tools_needed = ["samtools", "smalt", "spades.py"]
    for i in tools_needed:
        if not check_installed_tools(i, hard=False):
            print("Error! Can't run without install of {0} in PATH".format(i))
            sys.exit(1)
    logger = logging
    unittest.main()
