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
# I hate this line but it works :(
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))


from riboSeed.riboSim import ageSequence


sys.dont_write_bytecode = True

logger = logging


class RiboSimTest(unittest.TestCase):
    """
    """
    def setUp(self):
        self.fasta = os.path.join(
            os.path.dirname(__file__),
            "references", "riboSim_references", "original.fasta")
        self.temp_fasta = os.path.join(
            os.path.dirname(__file__),
            "references", "riboSim_references", "new.fasta")
        self.to_be_removed = []

    def test_Exes_bad_method(self):
        """check bad method arg given to Exes object"""
        random.seed(27)
        test_seq = str(
            "AAATAACAAAGACAAAAGACAAACAAAGAGAAGTAAACAATACATAAACAAAGACAAAAC" +
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
            "AAATAGAGAAACAATAAAATAAACAACACAACAAAAAGAAGAACACAAGAACAAAAAGAA")
        with open(self.fasta, "r") as infile:
            rec = list(SeqIO.parse(infile, "fasta"))[0]
        ageSequence(rec, freq=.1, end_length=60,
                    outfile=self.temp_fasta, logger=logger)
        with open(self.temp_fasta, "r") as infa:
            newrec = list(SeqIO.parse(infa, "fasta"))[0]
        self.assertEqual(test_seq, str(newrec.seq))

    def tearDown(self):
        """
        """
        for filename in self.to_be_removed:
            os.unlink(filename)
        pass


if __name__ == '__main__':
    unittest.main()
