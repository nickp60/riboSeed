# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas
The Goal of this is to have a unified place to put the useful
python 3.5 functions or templates

how I got the fastq file
# seqtk sample -s 27 ~/GitHub/FA/pseudochromosome/data/20150803_Abram1/ \
    reads/3123-1_1_trimmed.fastq .0005

bam file was from a riboseed mapping; md5: 939fbf2c282091aec0dfa278b05e94ec

mapped bam was made from bam file with the following command
 samtools view -Bh -F 4 /home/nicholas/GitHub/FB/Ecoli_comparative_genomics/
    scripts/riboSeed_pipeline/batch_coli_unpaired/map/
    mapping_20160906_region_7_riboSnag/
    test_smalt4_20160906_region_7_riboSnagS.bam >
     ~/GitHub/pyutilsnrw/tests/test_mapped.sam
md5: 27944249bf064ba54576be83053e82b0

"""

__version__ = "0.0.3"
import sys
import logging
import subprocess
import os
import unittest

from pyutilsnrw.utils3_5 import check_installed_tools, md5, file_len

from riboSeed.riboSeed2 import SeedGenome, ngsLib, LociCluster, LociMapping,\
    map_to_genome_smalt




sys.dont_write_bytecode = True

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class riboSeed2TestCase(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_riboseed2_tests")
        self.spades_dir = os.path.join(os.path.dirname(__file__),
                                       "output_riboseed_tests",
                                       "SPAdes_results")
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.ref_gb = os.path.join(self.ref_dir,
                                   'NC_011751.1.gb')
        self.ref_fasta = os.path.join(self.ref_dir,
                                      'cluster1.fasta')
        self.ref_Ffastq = os.path.join(self.ref_dir,
                                       'toy_reads1.fq')
        self.ref_Rfastq = os.path.join(self.ref_dir,
                                       'toy_reads2.fq')
        self.ref_bam_prefix = os.path.join(self.ref_dir,
                                           'test_bam_to_fastq')
        self.smalt_exe = "smalt"
        self.samtools_exe = "samtools"
        self.test_estimation_file = os.path.join(self.test_dir,
                                                 "est_distance.sam")
        self.map_results_prefix = os.path.join(self.test_dir,
                                               "test_mapping")
        self.fastq_results_prefix = os.path.join(self.test_dir,
                                                 "test_bam_to_fastq")
        self.to_be_removed = []

    def test_make_testing_dir(self):
        """ creates temp dir for all the files created in these tests
        """
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir)
        self.assertTrue(os.path.exists(self.test_dir))

    def test_LociMapping(self):
        testmapping = LociMapping(
            iteration=1,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"),
            finished=False)
        self.assertTrue(os.path.isdir(testmapping.mapping_subdir))

    def test_ngsLib(self):
        testlib = ngsLib(
            name="test",
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            readS0=None,
            ref_fasta=self.ref_fasta,
            smalt_dist_path=None,
            readlen=None,
            smalt_exe=self.smalt_exe)
        self.to_be_removed.append(testlib.smalt_dist_path)

    def test_SeedGenome(self):
        gen = SeedGenome(
            genbank_path=self.ref_gb,
            loci_clusters=None,
            output_root=self.test_dir)
        self.assertTrue(os.path.exists(os.path.join(self.test_dir,
                                                    "NC_011751.1.fasta")))
        print(gen.fasta_path)

    def test_map_to_genome_smalt(self):
        gen = SeedGenome(
            genbank_path=self.ref_gb,
            output_root=self.test_dir)
        gen.seq_ob = ngsLib(
            name="test",
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            readS0=None,
            ref_fasta=gen.fasta_path,
            smalt_dist_path=None,
            readlen=None,
            smalt_exe=self.smalt_exe)
        map_to_genome_smalt(
            seed_genome=gen,
            ngsLib=gen.seq_ob,
            map_results_prefix=gen.initial_map_prefix,
            cores=2,
            samtools_exe=self.samtools_exe,
            smalt_exe=self.smalt_exe,
            score_minimum=None,
            step=3, k=5,
            scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
            logger=logger)

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            os.unlink(filename)
        pass

if __name__ == '__main__':
    unittest.main()
    # tools_needed = ["samtools", "smalt", "spades.py"]
    # for i in tools_needed:
    #     if not check_installed_tools(i, hard=False):
    #         print("Error! Can't run without install of {0} in PATH".format(i))
    #         sys.exit(1)
    # logger = logging
