# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import shutil
# import subprocess
import os
import unittest
# import multiprocessing

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from argparse import Namespace
from unittest.mock import MagicMock, patch

# I hate this line but it works :(
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))


from riboSeed.shared_methods import md5

from riboSeed.riboSwap import remove_bad_contig, append_replacement_contigs
sys.dont_write_bytecode = True

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class riboSwapTestCase(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_riboSwap_tests")
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.ref_gb = os.path.join(self.ref_dir,
                                   'NC_011751.1.gb')
        self.bad_fasta = os.path.join(self.ref_dir,
                                      'test_de_fere_novo.fasta')
        self.good_fasta = os.path.join(self.ref_dir,
                                       'test_de_novo.fasta')
        self.ref_Ffastq = os.path.join(self.ref_dir,
                                       'toy_reads1.fq')
        self.ref_Rfastq = os.path.join(self.ref_dir,
                                       'toy_reads2.fq')
        self.ref_bam_prefix = os.path.join(self.ref_dir,
                                           'test_bam_to_fastq')
        self.smalt_exe = "smalt"
        self.bwa_exe = "bwa"
        self.samtools_exe = "samtools"
        self.spades_exe = "spades.py"
        self.quast_exe = "quast.py"
        self.python2_7_exe = "python2"
        self.test_estimation_file = os.path.join(self.test_dir,
                                                 "est_distance.sam")
        self.map_results_prefix = os.path.join(self.test_dir,
                                               "test_mapping")
        self.fastq_results_prefix = os.path.join(self.test_dir,
                                                 "test_bam_to_fastq")
        self.test_loci_file = os.path.join(os.path.dirname(__file__),
                                           str("references" + os.path.sep +
                                               'grouped_loci_reference.txt'))
        self.args = Namespace(skip_contol=False, kmers="21,33,55,77,99",
                              spades_exe="spades.py",
                              quast_exe="python2.7 quast.py",
                              cores=2)
        self.cores = 2
        self.maxDiff = 2000
        self.to_be_removed = []
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)
        # self.copy_fasta()

    def test_remove_bad_contig(self):
        outfile = os.path.join(self.test_dir, "contigs_minus_NODE_2.fasta")
        remove_bad_contig(infile=self.bad_fasta,
                          # outfile=self.good_fasta,
                          outfile=outfile,
                          # bad_name="NODE_3_length_222_cov_5.70526",
                          bad_name="NODE_3_",
                          logger=logger)
        with open(outfile, 'r') as of:
            recs = list(SeqIO.parse(of, 'fasta'))
            for i in recs:
                self.assertTrue("NODE_3_" not in i.id)
        self.to_be_removed.append(outfile)

    def test_fail_remove_bad_contig(self):
        outfile = os.path.join(self.test_dir, "contigs_minus_NODE_3.fasta")
        with self.assertRaises(ValueError):
            remove_bad_contig(infile=self.bad_fasta,
                              outfile=outfile,
                              bad_name="NODE_",
                              logger=logger)

    def test_append_replacement_contigs(self):
        outfile = os.path.join(self.test_dir, "contigs_minus_NODE_3.fasta")
        remove_bad_contig(infile=self.bad_fasta,
                          # outfile=self.good_fasta,
                          outfile=outfile,
                          # bad_name="NODE_3_length_222_cov_5.70526",
                          bad_name="NODE_3_",
                          logger=logger)
        append_replacement_contigs(infile=self.good_fasta, outfile=outfile,
                                   name_list="NODE_4_:NODE_5_".split(":"),
                                   logger=logger)
        self.assertEqual(md5(outfile),
                         md5(os.path.join(self.ref_dir,
                                          "ref_swapped_contigs.fasta")))
        self.to_be_removed.append(outfile)

    def test_fail_append_replacement_contigs(self):
        outfile = os.path.join(self.test_dir, "contigs_minus_NODE_3.fasta")
        remove_bad_contig(infile=self.bad_fasta,
                          outfile=outfile,
                          bad_name="NODE_3_",
                          logger=logger)
        with self.assertRaises(ValueError):
            append_replacement_contigs(
                infile=self.good_fasta, outfile=outfile,
                name_list="NODE_4_:NODE_5_:NODE_45_".split(":"),
                logger=logger)
        self.to_be_removed.append(outfile)

    def tearDown(self):
        """ delete temp files if no errors
        """
        if len(self.to_be_removed) != 0:
            for filename in self.to_be_removed:
                if os.path.exists(filename):
                    os.unlink(filename)

if __name__ == '__main__':
    unittest.main()
