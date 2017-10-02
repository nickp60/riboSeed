# -*- coding: utf-8 -*-
"""
@author: nicholas

"""
import sys
import shutil
import os
import unittest
import logging
import subprocess
from argparse import Namespace

from pyutilsnrw.utils3_5 import md5

from riboSeed.run_riboSeed import detect_or_create_config, parse_config, \
    new_log_for_diff

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 "with less than python 3.5")
class runRiboSeedTestCase(unittest.TestCase):
    """ Testing all the functions surrounding the actual plotting functions
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_run_riboSeed_tests")
        # self.ref_dir = os.path.join(
        #     os.path.dirname(__file__), "references", "")
        self.run_ref_dir = os.path.join(
            os.path.dirname(__file__),
            "references",
            "run_riboSeed_references", "")
        self.test_args = Namespace(
            outdir="",
            name="empty")

        os.makedirs(self.test_dir, exist_ok=True)
        self.to_be_removed = []

    def test_detect_or_create_config_exist(self):
        """ can we get the path to our config file back if it exists
        """
        path = detect_or_create_config(
            config_file=os.path.join(
                self.run_ref_dir, "sample_config.py"),
            theseargs=self.test_args,
            output_root="test", logger=logger)
        self.assertEqual(
            path, os.path.join(self.run_ref_dir, "sample_config.py"))

    def test_detect_or_create_config_noexist(self):
        """ can we get the path to our config file back if it dont exists
        """
        path = detect_or_create_config(
            config_file=os.path.join(
                self.run_ref_dir, "not_an_existing_sample_config.py"),
            output_root=self.test_dir,
            theseargs=self.test_args,
            newname="lookanewconfig",
            logger=logger)
        self.assertEqual(
            path, os.path.join(self.test_dir, "lookanewconfig.yaml"))

    def test_parse_config(self):
        """ can we read in values from our sample config file
        """
        conf = parse_config(
            config_file=os.path.join(
                self.run_ref_dir, "sample_config.py"),
            logger=logger)
        self.assertEqual(
            conf.MAUVE_ALIGNER_EXE,
            '/home/nicholas/mauve_snapshot_2015-02-13/linux-x64/mauveAligner')

    def test_new_log_for_diff(self):
        test_log = os.path.join(self.run_ref_dir, "run_riboSeed_notime.log")
        ref_timed_log = os.path.join(self.run_ref_dir, "sample_log.txt")
        ref_notime_log = os.path.join(self.run_ref_dir,
                                    "reference_run_riboSeed_notime.log")
        new_log_for_diff(ref_timed_log)
        self.assertEqual(
            md5(ref_notime_log), md5(test_log))
        self.to_be_removed.append(test_log)

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            try:
                os.unlink(filename)
            except IsADirectoryError:
                shutil.rmtree(filename)


if __name__ == '__main__':
    unittest.main()
