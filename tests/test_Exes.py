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

from riboSeed.classes import Exes


sys.dont_write_bytecode = True

logger = logging


@unittest.skipIf(shutil.which("bwa") is None or
                 shutil.which("quast") is None or
                 shutil.which("python2.7") is None or
                 shutil.which("spades.py") is None,
                 "bwa executable not found, skipping.If this isnt an " +
                 "error from travis deployment, you probably " +
                 "should install it")
class ExesTest(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.smalt_exe = "smalt"
        self.bwa_exe = "bwa"
        self.samtools_exe = "samtools"
        self.spades_exe = "spades.py"
        self.quast_exe = "quast"
        self.python_exe = "python"

    def test_Exes_bad_method(self):
        """check bad method arg given to Exes object"""
        with self.assertRaises(ValueError):
            Exes(samtools=self.samtools_exe,
                 quast=self.quast_exe,
                 smalt=self.smalt_exe,
                 python=self.python_exe,
                 spades=self.spades_exe,
                 bwa=self.bwa_exe,
                 method="bowtie")

    def test_Exes_bad_attribute(self):
        """ check bad instantiation of Exes object"""
        with self.assertRaises(AssertionError):
            Exes(samtools=None,
                 quast=self.quast_exe,
                 spades=self.spades_exe,
                 python=self.python_exe,
                 smalt=self.smalt_exe,
                 bwa=self.bwa_exe,
                 method="bwa")

    def test_Exes_bad_exe(self):
        """check Exes with nonexistant executable
        """
        with self.assertRaises(ValueError):
            Exes(samtools="nottheactualsamtools_exe",
                 quast=self.quast_exe,
                 python=self.python_exe,
                 smalt=self.smalt_exe,
                 spades=self.spades_exe,
                 bwa=self.bwa_exe,
                 method="bwa")

    def test_Exes_(self):
        """check with  executable"""
        test_exes = Exes(samtools=self.samtools_exe,
                         quast=self.quast_exe,
                         python=self.python_exe,
                         smalt=self.smalt_exe,
                         spades=self.spades_exe,
                         bwa=self.bwa_exe,
                         method="bwa")
        self.assertEqual(test_exes.mapper, shutil.which(test_exes.bwa))

    def tearDown(self):
        """
        """
        pass


if __name__ == '__main__':
    unittest.main()
