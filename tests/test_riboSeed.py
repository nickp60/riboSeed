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

from pyutilsnrw.utils3_5 import make_output_prefix, check_installed_tools,\
    copy_file, get_ave_read_len_from_fastq, get_number_mapped,\
    extract_mapped_and_mappedmates, keep_only_first_contig, md5,\
    combine_contigs, clean_temp_dir

from riboSeed.riboseed import  check_smalt_full_install,\
    map_to_ref_smalt, convert_bams_to_fastq



@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class utils3_5TestCase(unittest.TestCase):
    def setUp(self):
        self.ref_gb
        self.ref_fasta
        self.ref_Ffastq
        self.ref_Rfastq
        self.ref_disctance
        self.pileup
        pass

    def test_map_to_ref_smalt(self):
        map_to_ref_smalt(ref, ref_genome, fastq_read1, fastq_read2,
                         distance_results,
                         map_results_prefix, cores, samtools_exe,
                         smalt_exe, fastq_readS="",
                         read_len=100, step=3, k=5,
                         scoring="match=1,subst=-4,gapopen=-4,gapext=-3")

    def test_convert_bams_to_fastq(self):
        convert_bams_to_fastq(map_results_prefix,
                              fastq_results_prefix,
                              keep_unmapped)

    def test_run_spades(self):
        run_spades(output, ref, ref_as_contig, pe1_1='', pe1_2='', pe1_s='',
                   as_paired=True, keep_best=True, prelim=False,
                   groom_contigs='keep_first',
                   k="21,33,55,77,99", seqname='', spades_exe="spades.py")

    def test_check_samtools_pileup(self):
        check_samtools_pileup(self.pileup)

    def test_reconstruct_seq(self):
        reconstruct_seq(refpath, pileup, verbose=True, veryverb=False,
                    logger=None)

    def tearDown(self):
        pass

if __name__ == '__main__':
    samtools_exe = "samtools"
    check_smalt_full_install(smalt_exe, logger=None)
    unittest.main()
