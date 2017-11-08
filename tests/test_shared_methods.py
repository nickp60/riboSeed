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
# from riboSeed.shared_methods import  get_number_mapped



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
                                    "shared_method_references")
        self.test_md5s_prefix = os.path.join(self.ref_dir, "md5")
        self.test_bam_file = os.path.join(
            self.ref_dir, "mapping_reference_red.bam")
        self.startTime = time.time()
        self.cores = 2
        self.maxDiff = 2000
        self.to_be_removed = []

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            os.unlink(filename)
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))

if __name__ == '__main__':
    unittest.main()
