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

sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))

from riboSeed.run_riboSeed import detect_or_create_config, parse_config

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 "with less than python 3.5")
class riboSketchTestCase(unittest.TestCase):
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
        os.makedirs(self.test_dir, exist_ok=True)
        self.to_be_removed = []

    def test_detect_or_create_config_exist(self):
        """ can we get the path to our config file back if it exists
        """
        path = detect_or_create_config(
            config_file=os.path.join(
                self.run_ref_dir, "sample_config.py"),
            output_root="test", logger=logger)
        self.assertEqual(
            path, os.path.join(self.run_ref_dir, "sample_config.py"))

    def test_detect_or_create_config_noexist(self):
        """ can we get the path to our config file back if it exists
        """
        path = detect_or_create_config(
            config_file=os.path.join(
                self.run_ref_dir, "not_an_existing_sample_config.py"),
            output_root=self.test_dir,
            newname="lookanewconfig",
            logger=logger)
        self.assertEqual(
            path, os.path.join(self.test_dir, "lookanewconfig.py"))

    def test_parse_config(self):
        """ can we read in values from our sample config file
        """
        conf = parse_config(
            config_file=os.path.join(
                self.run_ref_dir, "sample_config.py"),
            logger=logger)
        self.assertEqual(
            conf.MAUVE_ALIGNER,
            '/home/nicholas/mauve_snapshot_2015-02-13/linux-x64/mauveAligner')

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
