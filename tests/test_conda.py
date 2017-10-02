# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import os
import unittest


logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class CondaInstallTests(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        pass

    def test_import_pandas(self):
        """ test pandas import
        """
        import pandas as pd

    def test_import_numpy(self):
        """ test pandas import
        """
        import numpy

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
