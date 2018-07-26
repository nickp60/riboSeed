# -*- coding: utf-8 -*-
"""
@author: nicholas

"""
import sys
import shutil
import os
import unittest
import time

from riboSeed.shared_methods import md5

from argparse import Namespace

sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))

from riboSeed.riboSketch import makeContigMoverCmds, findBestAlignments, \
    parseBackbones, parseDirContents, parseAlignmentDir, main, \
    plot_mauve_compare


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 "with less than python 3.5")
class riboSketchTestCase(unittest.TestCase):
    """ Testing all the functions surrounding the actual plotting functions
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_riboSketch_tests")
        self.ref_dir = os.path.join(
            os.path.dirname(__file__), "references", "")
        self.sketch_ref_dir = os.path.join(
            os.path.dirname(__file__),
            "references",
            "riboSketch_references", "")
        self.ref_png = os.path.join(
            os.path.dirname(__file__),
            "references",
            "riboSketch_references",
            "PrettyMauve.png")
        self.mauve_jar = os.path.join(
            "~", "mauve_snapshot_2015-02-13", "Mauve.jar")
        self.mauve_res_dir = os.path.join(
            os.path.dirname(__file__),
            "references",
            "riboSketch_references",
            "ref_vs_kleb_de_fere_novo", "alignment2"
            "")
        os.makedirs(self.test_dir, exist_ok=True)
        self.startTime = time.time()
        self.to_be_removed = []

    def test_parseDirContents(self):
        lst = parseDirContents(
            dirname=self.sketch_ref_dir, ref_ext="gb", assembly_ext="fasta")
        self.assertEqual(
            lst[0],
            os.path.join(self.sketch_ref_dir, "scannedScaffolds.gb"))
        self.assertEqual(
            lst[1][0],
            os.path.join(self.sketch_ref_dir, "mafft_msa.fasta"))

    def test_makeContigMoverCmds(self):
        cmds, results_path = makeContigMoverCmds(
            ref="test.gb", files=["fasta1.fa", "fasta2.fa"],
            outdir="reordering",
            mauve_jar="mauve.exe")
        cmd1 = "java -Xmx500m -cp mauve.exe org.gel.mauve.contigs.ContigOrderer -output reordering{0}ref_vs_fasta1 -ref test.gb -draft fasta1.fa".format(os.path.sep)
        cmd2 = "java -Xmx500m -cp mauve.exe org.gel.mauve.contigs.ContigOrderer -output reordering{0}ref_vs_fasta2 -ref test.gb -draft fasta2.fa".format(os.path.sep)

        self.assertEqual(cmds, [cmd1, cmd2])

    def test_findBestAlignments(self):
        testdir = os.path.join(self.sketch_ref_dir, "ref_vs_kleb_de_fere_novo")
        best_aln_dirs = findBestAlignments(testdir)
        self.assertEqual(
            os.path.join(testdir, "alignment2", ""),
            best_aln_dirs)

    def test_parseAlignmentDir(self):
        testdir = os.path.join(
            self.sketch_ref_dir, "ref_vs_kleb_de_fere_novo", "alignment2", "")
        assembly_list, backbone_list = parseAlignmentDir(dirlist=[testdir])
        self.assertEqual(assembly_list,
                         [os.path.join(testdir, "kleb_de_fere_novo.fasta")])
        self.assertEqual(backbone_list,
                         [os.path.join(testdir, "alignment2.backbone")])

    def test_parseBackbone(self):
        thisdir = os.path.join(
            self.sketch_ref_dir,
            "ref_vs_kleb_de_fere_novo", "alignment2", "alignment2.backbone")
        backbones = parseBackbones([thisdir])

        ref = [[25035, 34972, 3, 9940],
               [1, 21221, 9943, 31167],
               [36307, 36697, 31250, 31640],
               [36925, 50602, 31777, 45457],
               [50977, 52038, 45833, 46892],
               [52174, 65318, 47121, 60265],
               [70476, 80341, 60270, 70134],
               [85483, 95420, 70139, 80075],
               [100486, 105547, 80080, 85141],
               [21222, 25034, 0, 0],
               [34973, 36306, 0, 0],
               [36698, 36924, 0, 0],
               [50603, 50976, 0, 0],
               [52039, 52173, 0, 0],
               [65319, 70475, 0, 0],
               [80342, 85482, 0, 0],
               [95421, 100485, 0, 0],
               [0, 0, 31168, 31249],
               [0, 0, 31641, 31776],
               [0, 0, 45458, 45832],
               [0, 0, 46893, 47120]]
        self.assertEqual(backbones, [ref])

    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "Skipping this test on Travis CI. Too hard to debug")
    def test_plot_mauve_compare(self):
        backbone = os.path.join(
            self.mauve_res_dir, "alignment2.backbone")
        assembly_list = [
            os.path.join(self.mauve_res_dir, "kleb_de_fere_novo.fasta")]
        tempout = os.path.join(self.test_dir, "plot_out")
        os.makedirs(tempout, exist_ok=True)
        plot_mauve_compare(
            refgb=os.path.join(self.mauve_res_dir, "reference.gb"),
            assembly_list=assembly_list,
            backbones_list=[backbone],
            bufferlen=1000,
            breakwidth=100,
            aspect=.4,
            names=["reference", "kleb_de_fere_novo"],
            title="",
            output_prefix=os.path.join(tempout, "PrettyMauve"))
        self.assertEqual(md5(self.ref_png),
                         md5(os.path.join(tempout, "PrettyMauve.png")))
        self.to_be_removed.append(tempout)

    # @unittest.skipIf(not os.path.exists(
    #     os.path.join(
    #         os.path.dirname(
    #             os.path.dirname(
    #                 shutil.which("mauveAligner"))),
    #                  "Mauve.jar")),
    #                  "mauve jar not found, skipping." +
    #                  "If this isnt an error from travis, you " +
    #                  "probably should install it")
    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "Skipping this test on Travis CI.")
    def test_main(self):
        tempout = os.path.join(self.test_dir, "main_out")
        if  os.path.isdir(tempout):
            shutil.rmtree(tempout)
        self.args = Namespace(indir=self.mauve_res_dir,
                              outdir=tempout,
                              replot=False,
                              ref_ext="gb",
                              assembly_ext="fasta",
                              mauve_jar=self.mauve_jar,
                              names=None,
                              verbosity=2)
        main(self.args)
        self.assertEqual(md5(self.ref_png),
                         md5(os.path.join(tempout, "PrettyMauve.png")))
        self.to_be_removed.append(tempout)

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
