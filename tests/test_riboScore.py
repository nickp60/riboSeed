# -*- coding: utf-8 -*-
"""
@author: nicholas

"""
import sys
import logging
import shutil
import os
import unittest
import time

from unittest.mock import Mock


from riboSeed.shared_methods import md5

from riboSeed.riboScore import getSnagCmd, getSelectCmd, getScanCmd, \
    parseDirContents, make_nuc_nuc_recip_blast_cmds, merge_outfiles, \
    BLAST_tab_to_df, filter_recip_BLAST_df, checkBlastForMisjoin, \
    check_scan_select_snag_retruncodes

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
        self.score_ref_dir = os.path.join(
            os.path.dirname(__file__),
            "references",
            "riboScore_references", "")
        self.test_combine = os.path.join(
            os.path.dirname(__file__),
            "references",
            "riboScore_references",
            "test_combineA.tab")
        self.test_forward = os.path.join(
            os.path.dirname(__file__),
            "references",
            "riboScore_references",
            "forward.tab")
        self.test_reverse = os.path.join(
            os.path.dirname(__file__),
            "references",
            "riboScore_references",
            "reverse.tab")
        self.startTime = time.time()

        self.to_be_removed = []

    def test_parseDirContents(self):
        lst = parseDirContents(
            dirname=self.ref_dir, ref_ext="gb", assembly_ext="fasta")
        print(lst)

    def test_getScanCmd(self):
        res = getScanCmd(ref="test.fa", outroot="outdir", other_args="")
        res2 = getScanCmd(ref="test.gb", outroot="outdir", other_args="")
        ref_cmd = "ribo scan test.fa --min_length 5000 -o outdir{2}scan".format(
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
        res = getSelectCmd(gb="test.gb", outroot="outdir",
                           other_args="-s 16S:23S")
        ref_cmd = "ribo select test.gb -o outdir{2}select -s 16S:23S".format(
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
                         flank=20, outroot="outdir", other_args="")
        ref_cmd = "ribo snag test.gb clusters.txt -l 20 --just_extract -o outdir{2}snag".format(
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

    def test_merge_outfiles(self):
        """
        """
        merged_tab = merge_outfiles(
            filelist=[self.test_combine, self.test_combine],
            outfile=os.path.join(self.score_ref_dir, "temp_combined.tab"))
        self.assertEqual(
            md5(os.path.join(self.score_ref_dir, "test_combined.tab")),
            md5(merged_tab))
        self.to_be_removed.append(merged_tab)

    def test_single_merge_outfiles(self):
        """
        """
        merged_tab = merge_outfiles(
            filelist=[self.test_combine],
            outfile=os.path.join(self.score_ref_dir, "temp_combined.tab"))
        self.assertEqual(merged_tab, [self.test_combine])

    def test_BLAST_tab_to_df(self):
        colnames = [
            "query_id", "subject_id", "identity_perc", "alignment_length",
            "mismatches", "gap_opens", "q_start", "q_end", "s_start",
            "s_end", "evalue", "bit_score"]

        resultsdf = BLAST_tab_to_df(self.test_combine)
        self.assertEqual(resultsdf.columns.values.tolist(), colnames)

    def test_recip_blast(self):
        """ reciprocal blast testing.
        It doesnt really test much efficiently
        """
        df1 = BLAST_tab_to_df(self.test_forward)
        df2 = BLAST_tab_to_df(self.test_reverse)
        filtered_hits = filter_recip_BLAST_df(
            df1=df1,
            df2=df2,
            min_lens={"concatenated_genome_4001..10887": 500},
            min_percent=99.5,
            logger=logger)
        self.assertEqual(filtered_hits.shape, (2, 13))

    def test_checkBlastForMisjoin(self):
        df2 = BLAST_tab_to_df(self.test_reverse)
        flanking_hits = checkBlastForMisjoin(
            fasta="mock.fasta",
            df=df2,
            ref_lens={"concatenated_genome_4001..10887": 500},
            flanking=1000,
            BUF=50, logger=logger)
        self.assertEqual(
            flanking_hits[0],
            ["mock.fasta", "?",
             "NODE_1_length_105529_cov_19.8862_0_94652..101540_RC_",
             "concatenated_genome_4001..10887", "?"])

    def test_check_scan_select_snag_fail1(self):
        reslist = []
        for i in [1, 0, 0]:
            submock = Mock()
            submock.returncode = i
            reslist.append(submock)
        with self.assertRaises(SystemExit):
            check_scan_select_snag_retruncodes(
                subreturns=reslist, logger=logger)

    def test_check_scan_select_snag_fail2(self):
        reslist = []
        for i in [0, 1, 0]:
            submock = Mock()
            submock.returncode = i
            reslist.append(submock)
        with self.assertRaises(SystemExit):
            check_scan_select_snag_retruncodes(
                subreturns=reslist, logger=logger)

    def test_check_scan_select_snag_nofail(self):
        reslist = []
        for i in [0, 0, 1]:
            submock = Mock()
            submock.returncode = i
            reslist.append(submock)
        check_scan_select_snag_retruncodes(
            subreturns=reslist, logger=logger)

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            try:
                os.unlink(filename)
            except IsADirectoryError:
                shutil.rmtree(filename)
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))


if __name__ == '__main__':
    unittest.main()
