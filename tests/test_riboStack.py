# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import shutil
import os
import unittest

# from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import Namespace

# I hate this line but it works :(
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))


from pyutilsnrw.utils3_5 import md5
from riboSeed.riboStack import makeRegions, makeBedtoolsShuffleCmd, \
    samtoolsGetDepths

sys.dont_write_bytecode = True

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class riboStackTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_riboScore_tests")
        self.ref_dir = os.path.join(
            os.path.dirname(__file__), "references", "")
        self.stack_ref_dir = os.path.join(
            os.path.dirname(__file__),
            "references",
            "riboStack_references", "")
        self.test_combine = os.path.join(
            os.path.dirname(__file__),
            "references",
            "riboStack_references",
            "test_combineA.tab")
        self.gff = os.path.join(self.stack_ref_dir,
                                "scannedScaffolds.gff")

        self.to_be_removed = []

    def test_makeRegions_infer(self):
        rDNA_regions = os.path.join(self.stack_ref_dir, "rDNA_regions")

        makeRegions(outdir=self.stack_ref_dir, gff=self.gff, name="",
                    dest=rDNA_regions, logger=logger)
        # infer names
        with open(rDNA_regions, "r") as f:
            for i, line in enumerate(f):
                if i == 0:
                    self.assertEqual(line.split("\t")[0],
                                     "NC_000913.3_0")

    def test_makeRegions_name(self):
        gff = os.path.join(self.stack_ref_dir,
                           "scannedScaffolds.gff")
        rDNA_regions = os.path.join(self.stack_ref_dir, "rDNA_regions")

        makeRegions(outdir=self.stack_ref_dir, gff=gff, name="williewonka",
                    dest=rDNA_regions, logger=logger)
        # infer names
        with open(rDNA_regions, "r") as f:
            for i, line in enumerate(f):
                if i == 0:
                    self.assertEqual(line.split("\t")[0],
                                     "williewonka")
        self.to_be_removed.append(rDNA_regions)

    def test_makeBedtoolsShuffleCmd(self):
        ref_cmd = str(
            "bedtools shuffle -i rDNA_regions -g bedtools_genome > " +
            "shuffleRegions/sample_region_2")
        bedtools_cmds, bed_results_list = makeBedtoolsShuffleCmd(
            region="rDNA_regions",
            bedtools_exe="bedtools",
            destdir="shuffleRegions",
            n=2, genome="bedtools_genome")
        self.assertEqual(bedtools_cmds[1], ref_cmd)

    def test_samtoolsGetDepths(self):
        ref_cmd = str("samtools view -b -F 256 path/to/bam | samtools depth" +
                      " -b region_file - > sample_out_depth_1")

        ref_depth_path, sample_depths_paths, samtools_cmds = samtoolsGetDepths(
            samtools_exe="samtools",
            sample_file_list=["region_file"],
            bam="path/to/bam",
            ref_reg_file="region_file",
            outdir="" )
        self.assertEqual(samtools_cmds[1], ref_cmd)

    def tearDown(self):
        """ delete temp files
        """
        for filename in self.to_be_removed:
            try:
                os.unlink(filename)
            except IsADirectoryError:
                shutil.rmtree(filename)


if __name__ == '__main__':
    unittest.main()
