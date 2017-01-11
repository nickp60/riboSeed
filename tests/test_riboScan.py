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

# from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import Namespace
from unittest.mock import MagicMock, patch

# I hate this line but it works :(
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))


# from pyutilsnrw.utils3_5 import md5, file_len, copy_file, get_number_mapped
from riboSeed.riboScan import parse_fasta_header, make_barrnap_cmd, \
    add_locus_tags_to_gff

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
        self.ref_gb = os.path.join(self.ref_dir,
                                   'NC_011751.1.gb')
        self.ref_fasta = os.path.join(self.test_dir,
                                      'cluster1.fasta')
        self.good_contig = os.path.join(self.ref_dir,
                                        'contigs.fasta')
        self.short_contig = os.path.join(self.ref_dir,
                                         'contigs.fasta')
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
        # if not os.path.exists(self.test_dir):
        #     os.makedirs(self.test_dir, exist_ok=True)
        # self.copy_fasta()

    def test_gi_parse_fasta_header(self):
        """check with nonexistant executable"""
        header1 = str(">gi|218703261|ref|NC_011751.1| " +
                      "Escherichia coli UMN026 chromosome, complete genome\n")
        self.assertEqual("NC_011751.1", parse_fasta_header(header1))

    def test_nongi_parse_fasta_header(self):
        """check with nonexistant executable"""
        header2 = str(">testgenome gi|218703261|ref|NC_011751.1| " +
                      "Escherichia coli UMN026 chromosome, complete genome\n")
        self.assertEqual("testgenome", parse_fasta_header(header2))

    def test_fail_parse_fasta_header(self):
        """check with nonexistant executable"""
        header3 = str("testgenome|blabla|")
        with self.assertRaises(ValueError):
            parse_fasta_header(header3)

    @unittest.skipIf(shutil.which("barrnap") is None,
                     "barrnap executable not found. If this isnt an " +
                     "error from travis deployment, you probably " +
                     "should install it")
    def test_make_barrnap_cmd(self):
        cmd1 = make_barrnap_cmd(infasta="test.fasta", outgff="test.gff",
                                exe="barrnap", thresh=.2, kingdom='euk')
        ref_cmd1 = "{0} -kingdom euk test.fasta --reject 0.2 > test.gff".format(
            shutil.which("barrnap"))
        self.assertEqual(cmd1, ref_cmd1)

    def test_fail_exe_make_barrnap_cmd(self):
        with self.assertRaises(AssertionError):
            make_barrnap_cmd(infasta="test.fasta", outgff="test.gff",
                             exe="definitelynotbarrnap", thresh=.2,
                             kingdom='euk')

    @unittest.skipIf(shutil.which("barrnap") is None,
                     "barrnap executable not found. If this isnt an " +
                     "error from travis deployment, you probably " +
                     "should install it")
    def test_fail_thresh_make_barrnap_cmd(self):
        with self.assertRaises(AssertionError):
            make_barrnap_cmd(infasta="test.fasta", outgff="test.gff",
                             exe="barrnap", thresh=1.2,
                             kingdom='euk')

    # def test_add_locus_tags_to_gff(self):
    #     """check with  executable"""
    #     open_name = '%s.open' % __name__
    #     with patch(open_name, create=True) as mock_open:
    #         mock_open.return_value = MagicMock(spec=file)

    #     with open('/some/path', 'w') as f:
    #         f.write('something')
    #     file_handle = mock_open.return_value.__enter__.return_value
    #     file_handle.write.assert_called_with('something')

        # add_locus_tags_to_gff(gff, acc)

    # def tearDown(self):
    #     """ delete temp files if no errors
    #     """
    #     for filename in self.to_be_removed:
    #         os.unlink(filename)
    #     pass

if __name__ == '__main__':
    unittest.main()
