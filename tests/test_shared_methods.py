# -*- coding: utf-8 -*-
"""

"""
import sys
import logging
import shutil
import time
import os
import unittest

from .context import riboSeed
import riboSeed.shared_methods as sm
import riboSeed.classes as classes



@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class sharedMethodsTestCase(unittest.TestCase):
    """ tests for riboSeed's shared methods
    The tests that should be here have not actually been trasfered from the
    test modules for riboSnag yet.
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_shared_methods_tests")
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references",
                                    "shared_methods_references")
        self.test_fastq_file = os.path.join(self.ref_dir,
                                            'reads_reference.fastq')
        self.test_bam_file = os.path.join(self.ref_dir,
                                          "mapping_reference.bam")
        self.test_md5s_prefix = os.path.join(self.ref_dir, "md5")
        self.test_bam_file = os.path.join(
            self.ref_dir, "mapping_reference.bam")
        self.test_multifasta = os.path.join(
            self.ref_dir,
            "test_multiseqs_reference.fasta")
        self.test_singlefasta = os.path.join(
            self.ref_dir,
            "test_only_first_reference.fasta")
        self.genbank_filename = os.path.join(
            self.ref_dir,
            'single.gb')
        self.multigenbank_filename = os.path.join(
            self.ref_dir,
            'multi.gb')
        self.test_md5s_prefix = os.path.join(
            self.ref_dir, "md5")
        self.startTime = time.time()
        self.cores = 2
        self.maxDiff = 2000
        self.to_be_removed = []

    def test_combine_contigs(self):
        pass

    def test_average_read_len(self):
        """ tests get_ave_read_len_from_fastq
        this probably could/should be refined to have a better test
        """
        mean_read_len = classes.get_ave_read_len_from_fastq(self.test_fastq_file, N=5)
        self.assertEqual(217.8, mean_read_len)

    @unittest.skipIf(shutil.which("samtools") is None,
                     "samtools executable not found, skipping." +
                     "If this isnt an error from travis deployment, you " +
                     "probably should install it")
    def test_get_number_mapped(self):
        """ checks flagstat
        """
        result = sm.get_number_mapped(self.test_bam_file, shutil.which("samtools"))
        reference = "151 + 0 mapped (0.56% : N/A)"
        self.assertEqual(result, reference)

    def test_keep_only_first_contig(self):
        pass

    def test_get_fasta_lengths(self):
        """ get the lengths of the multifasta entries
        """
        self.assertEqual(sm.get_fasta_lengths(self.test_singlefasta), [169])
        self.assertEqual(sm.get_fasta_lengths(self.test_multifasta),
                         [169, 161, 159, 159, 151, 133, 128])
    def test_file_len(self):
        """ test against file of known length
        """
        self.assertEqual(sm.file_len(self.test_singlefasta), 4)

    def test_multisplit(self):
        """ split a string that has multiple delimiters
        """
        test_string = "look_this+is+a locus_that_is+multi-delimited"
        list_of_things = sm.multisplit(["-", "_", "+", " "], test_string)
        test_other_string = "look_this+is+a\faillocus_that_is+multi-delimited"
        list_of_other_things = sm.multisplit(["-", "_", "+", " "],
                                          test_other_string)
        self.assertEqual(list_of_things, ["look", "this", "is", "a", "locus",
                                          "that", "is", "multi", "delimited"])
        self.assertNotEqual(list_of_other_things, ["look", "this", "is", "a",
                                                   "locus", "that", "is",
                                                   "multi", "delimited"])

    def test_get_genbank_record(self):
        """Reads records from a GenBank file.
        """
        records = sm.get_genbank_record(self.genbank_filename)
        assert isinstance(records, list)
        multirecords = sm.get_genbank_record(self.multigenbank_filename,
                                          first_only=False)
        assert isinstance(multirecords, list)

    def test_md5_strings(self):
        """ minimal md5 examples with strings
        """
        self.assertEqual(sm.md5("thisstringisidenticalto", string=True),
                         sm.md5("thisstringisidenticalto", string=True))
        self.assertNotEqual(sm.md5("thisstringisntidenticalto", string=True),
                            sm.md5("thisstringisnotidenticalto", string=True))

    def test_md5_files(self):
        """ test file contests identiy
        """
        md5_a = sm.md5(str(self.test_md5s_prefix + "_a.txt"))
        md5_b = sm.md5(str(self.test_md5s_prefix + "_b.txt"))
        md5_fail = sm.md5(str(self.test_md5s_prefix + "_fail.txt"))
        self.assertEqual(md5_a, md5_b)
        self.assertNotEqual(md5_a, md5_fail)

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            os.unlink(filename)
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))

if __name__ == '__main__':
    unittest.main()
