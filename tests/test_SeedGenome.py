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
from Bio import SeqIO

from .context import riboSeed
from riboSeed.riboSeed import SeedGenome, add_coords_to_clusters
from riboSeed.shared_methods import parse_clustered_loci_file


logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class SeedGenomeTest(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_SeedGenome_tests")
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.ref_gb = os.path.join(self.ref_dir,
                                   'NC_011751.1.gb')
        self.ref_tiny_gb = os.path.join(self.ref_dir,
                                        'scannedScaffolds.gb')
        self.ref_fasta = os.path.join(self.test_dir,
                                      'cluster1.fasta')
        self.test_loci_file = os.path.join(os.path.dirname(__file__),
                                           str("references" + os.path.sep +
                                               'grouped_loci_reference.txt'))
        self.to_be_removed = []
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)
        self.copy_fasta()

    def copy_fasta(self):
        """ make a disposable copy of the cluster fasta
        """
        shutil.copy(os.path.join(self.ref_dir, 'cluster1.fasta'),
                    self.ref_fasta)
        self.to_be_removed.append(self.ref_fasta)

    def test_SeedGenome(self):
        """ Can we create a seedGenome object
        """
        with open(self.ref_tiny_gb, "r") as gbf:
            rec = next(SeqIO.parse(gbf, "genbank"))
        gen = SeedGenome(
            max_iterations=2,
            clustered_loci_txt=self.test_loci_file,
            genbank_path=self.ref_tiny_gb,
            loci_clusters=None,
            output_root=self.test_dir)
        self.assertTrue(os.path.exists(os.path.join(self.test_dir,
                                                    "scannedScaffolds.fasta")))
        self.assertTrue(os.path.exists(self.test_dir))
        self.assertEqual(tuple(gen.seq_records)[0].id, rec.id)

    def test_SeedGenome_bad_instatiation(self):
        """ Does SeedGenome fail when missing attributes """
        with self.assertRaises(ValueError):
            SeedGenome(
                max_iterations=2,
                genbank_path=self.ref_gb,
                clustered_loci_txt=self.test_loci_file,
                loci_clusters=None,
                output_root=None)

    def test_add_coords_to_SeedGenome(self):
        """ can we parse a loci_coords file and add to seedGenome
        BTW, this reveals some pokiness in the parse_clustered_loci_file
        method (~20 seconds), and in add_coords_to_clusters (~6 seconds).
        Pull requests welcome :)

        also testing pad genbank to save a bit of time, as setup is the same
        """
        gen = SeedGenome(
            max_iterations=1,
            genbank_path=self.ref_gb,
            clustered_loci_txt=self.test_loci_file,
            output_root=self.test_dir,
            logger=logger)
        gen.loci_clusters = parse_clustered_loci_file(
            filepath=gen.clustered_loci_txt,
            gb_filepath=gen.genbank_path,
            output_root=self.test_dir,
            padding=100,
            circular=False,
            logger=logger)
        add_coords_to_clusters(seedGenome=gen, logger=logger)
        self.assertEqual(
            gen.loci_clusters[0].loci_list[0].start_coord, 4656045)
        self.assertEqual(
            gen.loci_clusters[0].loci_list[0].end_coord, 4657586)
        gen.pad_genbank(pad=1000, circular=True, logger=logger)
        with open(self.ref_gb, "r") as ingb:
            oldrec = list(SeqIO.parse(ingb, "genbank"))[0]
        with open(gen.ref_fasta, "r") as infa:
            newrec = list(SeqIO.parse(infa, "fasta"))[0]
        self.assertEqual(str(oldrec.seq),
                         str(newrec.seq[1000: - 1000]))

    def test_purge_old_files(self):
        gen = SeedGenome(
            max_iterations=1,
            genbank_path=self.ref_gb,
            clustered_loci_txt=self.test_loci_file,
            output_root=self.test_dir,
            logger=logger)
        target = os.path.join(self.ref_dir, 'thisbettergetnuked.fasta')
        shutil.copy(os.path.join(self.ref_dir, 'cluster1.fasta'), target)
        self.assertTrue(os.path.isfile(target))
        gen.iter_mapping_list[0].s_map_bam = target
        gen.purge_old_files(all_iters=True, logger=logger)
        self.assertFalse(os.path.isfile(target))

    def test_early_purge(self):
        gen = SeedGenome(
            max_iterations=1,
            genbank_path=self.ref_gb,
            clustered_loci_txt=self.test_loci_file,
            output_root=self.test_dir,
            logger=logger)
        gen.this_iteration = 1
        with self.assertRaises(AssertionError):
            gen.purge_old_files(all_iters=False, logger=logger)

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            os.unlink(filename)
        pass


if __name__ == '__main__':
    unittest.main()
