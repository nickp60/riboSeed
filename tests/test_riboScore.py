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


from pyutilsnrw.utils3_5 import md5
from riboSeed.riboScore import getSnagCmd, getSelectCmd, getScanCmd, \
    parseDirContents, make_nuc_nuc_recip_blast_cmds
sys.dont_write_bytecode = True

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class riboScoreTestCase(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_riboScore_tests")
        self.ref_dir = os.path.join(
            os.path.dirname(__file__), "references", "")
        self.scan_ref_dir = os.path.join(os.path.dirname(__file__),
                                         "references",
                                         "riboScan_references")
        self.to_be_removed = []

    def test_parseDirContents(self):
        lst = parseDirContents(
            dirname=self.ref_dir, ref_ext="gb", assembly_ext="fasta")
        print(lst)

    def test_getScanCmd(self):
        res = getScanCmd(ref="test.fa", outroot="outdir")
        res2 = getScanCmd(ref="test.gb", outroot="outdir")
        ref_cmd = "{0} {1} test.fa --min_length 5000 -o outdir{2}scan".format(
            sys.executable,
            os.path.join("..", "..",
                         os.path.dirname(os.path.dirname(__file__)),
                         "riboSeed",
                         "riboScan.py"),
            os.path.sep)
        self.assertEqual(
            res[0],
            ref_cmd)
        self.assertEqual(
            res[1], os.path.join("outdir", "scan", "scannedScaffolds.gb"))
        self.assertEqual(res2[0], None)

    def test_getSelectCmd(self):
        res = getSelectCmd(gb="test.gb", outroot="outdir")
        ref_cmd = "{0} {1} test.gb -s 16S:23S -o outdir{2}select".format(
            sys.executable,
            os.path.join(
                "..", "..",
                os.path.dirname(os.path.dirname(__file__)),
                "riboSeed",
                "riboSelect.py"),
            os.path.sep)
        self.assertEqual(res[0], ref_cmd)
        self.assertEqual(
            res[1],
            os.path.join("outdir", "select", "riboSelect_grouped_loci.txt"))

    def test_getSnagCmd(self):
        res = getSnagCmd(scangb="test.gb", cluster="clusters.txt",
                         flank=20, outroot="outdir")
        ref_cmd = "{0} {1} test.gb clusters.txt -l 20 --just_extract -o outdir{2}snag".format(
            sys.executable,
            os.path.join(
                "..", "..",
                os.path.dirname(os.path.dirname(__file__)),
                "riboSeed",
                "riboSnag.py"),
            os.path.sep)
        self.assertEqual(res[0], ref_cmd)

    def test_make_nuc_nuc_recip_blast_cmds(self):
        cmds, forward, recip = make_nuc_nuc_recip_blast_cmds(
            query_list=["assembly1.fasta", "assembly2.fasta"],
            output="outdir", subject_file="reference", logger=logger)
        self.assertEqual(
            cmds[0],
            "blastn -out outdir/assembly1_vs_ref.tab -outfmt 6 " +
            "-query assembly1.fasta -subject reference -num_threads 1 " +
            "-num_alignments 50")

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            try:
                os.unlink(filename)
            except IsADirectoryError:
                shutil.rmtree(filename)
        pass

if __name__ == '__main__':
    unittest.main()
