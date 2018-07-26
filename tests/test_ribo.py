# -*- coding: utf-8 -*-
"""
@author: nicholas

"""
import sys
import logging
import shutil
import os
import unittest
import time

from unittest.mock import Mock


from riboSeed.ribo import main

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class riboTestCase(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.startTime = time.time()

    def test_main_return_1_invalid(self):
        args = ["ribo", "notavalidarg"]
        self.assertEqual(main(args), 1)

    def test_main_return_0_version(self):
        args = ["ribo", "-v"]
        self.assertEqual(main(args), 0)

    def test_main_return_0_help(self):
        args = ["ribo"]
        self.assertEqual(main(args), 0)

    def tearDown(self):
        """
        """
        pass


if __name__ == '__main__':
    unittest.main()
