# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import os
import unittest


from riboSeed.riboSeed import NgsLib, nonify_empty_lib_files

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class NgsLibTest(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_NgsLib_tests")
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.ref_fasta = os.path.join(self.test_dir,
                                      'cluster1.fasta')
        self.ref_Ffastq = os.path.join(self.ref_dir,
                                       'toy_reads1.fq')
        self.ref_Rfastq = os.path.join(self.ref_dir,
                                       'toy_reads2.fq')
        self.smalt_exe = "smalt"
        self.bwa_exe = "bwa"
        self.to_be_removed = []
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)

    def test_NgsLib(self):
        """ Can we create an NgsLib object correctly
        """
        # make a non-master object
        testlib_pe_s = NgsLib(
            name="test",
            master=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            readS0="dummy",
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)
        testlib_s = NgsLib(
            name="test",
            master=True,
            readF=None,
            readR=None,
            readS0=self.ref_Ffastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)
        self.assertEqual(testlib_s.libtype, "s_1")
        self.assertEqual(testlib_s.readlen, 145.0)
        self.assertEqual(testlib_s.liblist, [self.ref_Ffastq])
        # test unnamed fails
        with self.assertRaises(ValueError):
            NgsLib(
                name=None,
                master=False,
                readF=self.ref_Ffastq,
                readR=self.ref_Rfastq,
                readS0="dummy",
                ref_fasta=self.ref_fasta,
                mapper_exe=self.smalt_exe)
        self.assertEqual(testlib_pe_s.libtype, "pe_s")
        self.assertEqual(testlib_pe_s.readlen, None)
        # test fails with singe PE file
        with self.assertRaises(ValueError):
            NgsLib(
                name=None,
                master=False,
                readF=self.ref_Ffastq,
                readR=None,
                readS0="dummy",
                ref_fasta=self.ref_fasta,
                mapper_exe=self.smalt_exe)

        # check master files cannot bge deleted
        testlib_pe = NgsLib(
            name="test",
            master=True,
            make_dist=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe,
            logger=logger)
        self.assertEqual(
            1,  # return code for no deleting tool place
            testlib_pe.purge_old_files(master=testlib_pe_s, logger=logger))
        self.assertTrue(os.path.isfile(self.ref_Ffastq))
        self.assertTrue(os.path.isfile(self.ref_Rfastq))

        # test killer lib that tries to purge files that are in a master ob
        testlib_killer = NgsLib(
            name="test",
            master=False,
            make_dist=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe,
            logger=logger)
        self.assertEqual(
            1,  # return code for no deleting tool place
            testlib_killer.purge_old_files(master=testlib_pe_s, logger=logger))
        self.assertTrue(os.path.isfile(self.ref_Ffastq))
        self.assertTrue(os.path.isfile(self.ref_Rfastq))

    def test_dont_check_nonmaster_read_len(self):
        testlib_pe = NgsLib(
            name="test",
            master=False,
            make_dist=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe,
            logger=logger)
        self.assertEqual(testlib_pe.readlen, None)

    def test_lib_check(self):
        """ does the NgsLib identify empty libraries
        """
        empty_file = os.path.join(self.test_dir, "test_not_real_file")
        # make an empty file
        with open(empty_file, 'w') as ef:
            pass
        ngs_ob = NgsLib(
            name="test",
            master=False,
            readF=self.ref_Ffastq,
            readR=empty_file,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)
        nonify_empty_lib_files(ngsLib=ngs_ob, logger=logger)
        self.assertTrue(ngs_ob.readR is None)
        self.to_be_removed.append(empty_file)

    def tearDown(self):
        """
        """
        for filename in self.to_be_removed:
            os.unlink(filename)
        pass
        pass


if __name__ == '__main__':
    unittest.main()
