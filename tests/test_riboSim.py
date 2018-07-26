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
import random
from Bio import SeqIO
from argparse import Namespace
# I hate this line but it works :(
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))


from riboSeed.riboSim import ageSequence, main
from riboSeed.shared_methods import md5



sys.dont_write_bytecode = True

logger = logging


class RiboSimTest(unittest.TestCase):
    """
    """
    def setUp(self):
        self.results_dir = os.path.join(os.path.dirname(__file__),
                                        "output_riboSim_tests")
        self.fasta = os.path.join(
            os.path.dirname(__file__),
            "references", "riboSim_references", "original.fasta")
        self.temp_fasta = os.path.join(
            os.path.dirname(__file__),
            "references", "riboSim_references", "new.fasta")
        self.temp_fasta2 = os.path.join(
            os.path.dirname(__file__),
            "references", "riboSim_references", "new2.fasta")
        self.mut_fasta = os.path.join(
            os.path.dirname(__file__),
            "references", "riboSim_references", "mutated_ten_seed_27.fasta")
        self.to_be_removed = []

    def test_correct_mutations(self):
        """check the aged sequence"""
        test_seq = str(
            "AAAAAAAAATAACAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAT" +
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
            "AAAAAGTAAAAAAAAAAAAAAAAAAATAAAGAAAGAAAAAAAAAAAAAAAAAAAACAAGA")
        with open(self.fasta, "r") as infile:
            rec = list(SeqIO.parse(infile, "fasta"))[0]
        ageSequence(rec, freq=.1, end_length=60,
                    outfile=self.temp_fasta, seed=27, logger=logger)
        with open(self.temp_fasta, "r") as infa:
            newrec = list(SeqIO.parse(infa, "fasta"))[0]
        self.assertEqual(test_seq, str(newrec.seq))
        self.to_be_removed.append(self.temp_fasta)

    def test_correct_mutation_rate(self):
        """check number of substitutions in aged sequence is sensible"""
        with open(self.fasta, "r") as infile:
            rec = list(SeqIO.parse(infile, "fasta"))[0]
        ageSequence(rec, freq=.1, end_length=0,
                    outfile=self.temp_fasta2, seed=27, logger=logger)
        with open(self.temp_fasta2, "r") as infa:
            newrec = list(SeqIO.parse(infa, "fasta"))[0]
        self.assertEqual(newrec.seq.count("A"), 162)
        self.to_be_removed.append(self.temp_fasta2)

    def test_main_integration(self):
        if os.path.exists(self.results_dir):
            shutil.rmtree(self.results_dir)
        args = Namespace(output=self.results_dir,
                         fasta=self.fasta,
                         verbosity=2, frequency=.1,
                         end_length=0, seed=27)
        main(args)
        self.assertEqual(
            md5(self.mut_fasta),
            md5(os.path.join(self.results_dir, "original.fasta")))
        self.to_be_removed.append(
            os.path.join(self.results_dir, "original.fasta"))

    def tearDown(self):
        """
        """
        for filename in self.to_be_removed:
            try:
                os.unlink(filename)
            except IsADirectoryError:
                shutil.rmtree(filename)


if __name__ == '__main__':
    unittest.main()
