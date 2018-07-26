# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import shutil
import time
# import subprocess
import os
import unittest
# import multiprocessing

# from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import Namespace
from unittest.mock import MagicMock, patch

from riboSeed.shared_methods import md5
from riboSeed.riboScan import parse_fasta_header, \
    add_locus_tags_to_gff, combine_gbs, append_accession_and_version, \
    make_seqret_cmd, splitMultifasta, getFastas, checkSingleFasta, main
from riboSeed.shared_methods import make_barrnap_cmd

sys.dont_write_bytecode = True

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class riboSeedTestCase(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_riboScan_tests")
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.scan_ref_dir = os.path.join(os.path.dirname(__file__),
                                         "references",
                                         "riboScan_references")
        self.no_locus_gff = os.path.join(self.scan_ref_dir,
                                         "no_locus.gff")
        self.with_locus_gff = os.path.join(self.scan_ref_dir,
                                           "with_locus.gff")
        self.combined_file = os.path.join(self.scan_ref_dir,
                                          "combined_gbs_just_kidding.txt")
        self.multifasta = os.path.join(self.scan_ref_dir,
                                       "multiFasta.fasta")
        self.no_accession_gb = os.path.join(self.scan_ref_dir,
                                            'no_accession_or_version.gb')
        self.with_accession_gb = os.path.join(self.scan_ref_dir,
                                              'with_accession_or_version.gb')

        self.good_contig = os.path.join(self.ref_dir,
                                        'contigs.fasta')
        self.short_contig = os.path.join(self.ref_dir,
                                         'contigs.fasta')
        self.args = Namespace(skip_contol=False, kmers="21,33,55,77,99",
                              spades_exe="spades.py",
                              quast_exe="python2.7 quast.py",
                              cores=2)
        self.startTime = time.time()
        self.cores = 2
        self.maxDiff = 2000
        self.to_be_removed = []

    def test_gi_parse_fasta_header(self):
        """check valid fasta header parsing"""
        header1 = str(">gi|218703261|ref|NC_011751.1| " +
                      "Escherichia coli UMN026 chromosome, complete genome\n")
        self.assertEqual("NC_011751.1", parse_fasta_header(header1))

    def test_nongi_parse_fasta_header(self):
        """check nonstandard fasta header parsing"""
        header2 = str(">testgenome gi|218703261|ref|NC_011751.1| " +
                      "Escherichia coli UMN026 chromosome, complete genome\n")
        self.assertEqual("testgenome", parse_fasta_header(header2))

    def test_fail_parse_fasta_header(self):
        """ensudre bad fastas fail header parsing"""
        header3 = str("testgenome|blabla|")
        with self.assertRaises(ValueError):
            parse_fasta_header(header3)

    def test_make_barrnap_cmd(self):
        """ check barrnap command construction  """
        cmd1 = make_barrnap_cmd(infasta="test.fasta", outgff="test.gff",
                                exe="barrnap.exe", thresh=.2, kingdom='euk')
        ref_cmd1 = str("barrnap.exe -k euk test.fasta --reject 0.2 --threads" +
                       " 1 --evalue 1e-06 > test.gff")
        self.assertEqual(cmd1, ref_cmd1)

    def test_fail_barrnap_exe_main(self):
        """ fail main with bad barrnap exe"""
        with self.assertRaises(AssertionError):
            main(Namespace(contigs="test.fasta",
                           output="test",
                           barrnap_exe="definitelynotbarrnap",
                           cores=2))

    # @unittest.skipIf(shutil.which("barrnap") is None,
    #                  "barrnap executable not found. If this isnt an " +
    #                  "error from travis deployment, you probably " +
    #                  "should install it")
    def test_fail_thresh_make_barrnap_cmd(self):
        with self.assertRaises(AssertionError):
            make_barrnap_cmd(infasta="test.fasta", outgff="test.gff",
                             exe="barrnap", thresh=1.2,
                             kingdom='euk')

    def test_add_locus_tags_to_gff(self):
        """check create proper gff """
        add_locus_tags_to_gff(gff=self.no_locus_gff,
                              acc="BA000007.2")
        new_gff = str(os.path.splitext(self.no_locus_gff)[0] + "_tagged.gff")
        self.assertEqual(md5(new_gff), md5(self.with_locus_gff))
        self.to_be_removed.append(new_gff)

    def test_combine_gbs(self):
        """ make sure we combine genbank files properly"""
        temp_gb = os.path.join(self.scan_ref_dir, "temp_combined.gb")
        combine_gbs(finalgb=temp_gb,
                    gb_list=[self.no_locus_gff,
                             self.with_locus_gff])
        self.assertEqual(md5(temp_gb), md5(self.combined_file))
        # self.to_be_removed.append(temp_gb)

    def test_append_accession_and_version(self):
        temp_gb2 = os.path.join(self.scan_ref_dir, "temp_acessioned.gb")
        append_accession_and_version(accession="BA000007.2",
                                     ingb=self.no_accession_gb,
                                     finalgb=temp_gb2)
        self.assertEqual(md5(temp_gb2), md5(self.with_accession_gb))
        self.to_be_removed.append(temp_gb2)

    @unittest.skipIf(
        shutil.which("seqret") is None,
        "seqret executable not found, skipping." +
        "If this isnt an error from travis deployment, you probably " +
        "should install it")
    def test_make_seqret_cmd(self):
        test_cmd = make_seqret_cmd(exe="seqret",
                                   outgb="dest_file.gb",
                                   infasta="input_sequence.fasta",
                                   ingff="input_annos.gff")
        ref_cmd = str(
            "{0} -sequence input_sequence.fasta -feature -fformat gff3 " +
            "-fopenfile input_annos.gff -osformat genbank -auto " +
            "-outseq dest_file.gb").format(shutil.which("seqret"))
        self.assertEqual(ref_cmd, test_cmd)

    def test_getFastas(self):
        getFastas(inp=self.multifasta, output_root=self.test_dir,
                  name="snorkel", logger=logger)
        self.assertTrue(os.path.isfile(
            os.path.join(self.test_dir, "contigs", "snorkel_1.fa")))
        with self.assertRaises(SystemExit):
            getFastas(inp="notanactualfile", output_root=self.test_dir,
                      name="snorkel", logger=logger)
        self.to_be_removed.append(os.path.join(self.test_dir, "contigs"))

    def test_checkSingleFasta(self):
        """  ensure failure if a multifasta snuck through"""
        with self.assertRaises(SystemExit):
            checkSingleFasta(self.multifasta, logger=logger)

    def tearDown(self):
        """ delete temp files if no errors """
        for filename in self.to_be_removed:
            try:
                os.unlink(filename)
            except IsADirectoryError:
                shutil.rmtree(filename)
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))

if __name__ == '__main__':
    unittest.main()
