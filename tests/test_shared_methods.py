# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import shutil
import time
import os
import unittest


from riboSeed.shared_methods import multisplit, md5, get_number_mapped



@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class riboSeedTestCase(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_shared_methods_tests")
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references",
                                    "shared_method_references")
        self.test_md5s_prefix = os.path.join(self.ref_dir, "md5"))
        self.test_bam_file = os.path.join(
            self.ref_dir, "mapping_reference_red.bam")
        self.startTime = time.time()
        self.cores = 2
        self.maxDiff = 2000
        self.to_be_removed = []

    def test_multisplit(self):
        """ split a string that has multiple delimiters
        """
        test_string = "look_this+is+a locus_that_is+multi-delimited"
        list_of_things = multisplit(["-", "_", "+", " "], test_string)
        test_other_string = "look_this+is+a\faillocus_that_is+multi-delimited"
        list_of_other_things = multisplit(["-", "_", "+", " "],
                                          test_other_string)
        self.assertEqual(list_of_things, ["look", "this", "is", "a", "locus",
                                          "that", "is", "multi", "delimited"])
        self.assertNotEqual(list_of_other_things, ["look", "this", "is", "a",
                                                   "locus", "that", "is",
                                                   "multi", "delimited"])
    def test_md5_strings(self):
        """ minimal md5 examples with strings
        """
        self.assertEqual(md5("thisstringisidenticalto", string=True),
                         md5("thisstringisidenticalto", string=True))
        self.assertNotEqual(md5("thisstringisntidenticalto", string=True),
                            md5("thisstringisnotidenticalto", string=True))

    def test_md5_files(self):
        """ test file contests identiy
        """
        md5_a = md5(str(self.test_md5s_prefix + "_a.txt"))
        md5_b = md5(str(self.test_md5s_prefix + "_b.txt"))
        md5_fail = md5(str(self.test_md5s_prefix + "_fail.txt"))
        self.assertEqual(md5_a, md5_b)
        self.assertNotEqual(md5_a, md5_fail)

    @unittest.skipIf(shutil.which("samtools") is None,
                     "samtools executable not found, skipping." +
                     "If this isnt an error from travis deployment, you " +
                     "probably should install it")
    def test_get_number_mapped(self):
        """ checks flagstat
        """
        result = get_number_mapped(self.test_bam_file, self.samtools_exe)
        reference = "2 + 0 mapped (0.67% : N/A)"
        self.assertEqual(result, reference)
