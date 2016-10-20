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
import time
import sys
import shutil
import logging
import subprocess
import os
import unittest
import hashlib
import glob
import argparse
from Bio import SeqIO
sys.dont_write_bytecode = True


from pyutilsnrw.utils3_5 import get_genbank_seq, get_genbank_record

from riboSeed.riboSnag import parse_clustered_loci_file, \
    extract_coords_from_locus, strictly_increasing, \
    stitch_together_target_regions, get_genbank_rec_from_multigb,\
    pad_genbank_sequence


logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class riboSnag_TestCase(unittest.TestCase):
    def setUp(self):
        self.curdir = os.getcwd()
        self.testdirname = os.path.join(os.path.dirname(__file__),
                                        "output_utils3_5_tests")
        self.test_loci_file = os.path.join(os.path.dirname(__file__),
                                           str("references" + os.path.sep +
                                               'grouped_loci_reference.txt'))
        self.test_gb_file = os.path.join(os.path.dirname(__file__),
                                         str("references" + os.path.sep +
                                             'NC_011751.1.gb'))
        self.test_loci_file = os.path.join(os.path.dirname(__file__),
                                           str("references" + os.path.sep +
                                               'grouped_loci_reference.txt'))
        self.test_cluster1 = os.path.join(os.path.dirname(__file__),
                                          str("references" + os.path.sep +
                                              'cluster1.fasta'))
        self.samtools_exe = "samtools"

    def test_parse_loci(self):
        """
        """
        clusters = parse_clustered_loci_file(self.test_loci_file,
                                             logger=logger)
        ref_loci_list = ['NC_011751.1',
                         ['ECUMN_16S_6', 'ECUMN_23S_6', 'ECUMN_5S_7']]

        self.assertEqual(clusters[0], ref_loci_list)

    def test_get_genbank_seq_matching_id(self):
        records = get_genbank_record(self.test_gb_file)
        record = get_genbank_rec_from_multigb(recordID='NC_011751.1',
                                              genbank_records=records)
        self.assertEqual(records[0].seq, record.seq)

    def test_extract_coords_from_locus(self):
        records = get_genbank_record(self.test_gb_file)
        coord_list = extract_coords_from_locus(record=records[0],
                                               locus_tag_list=["ECUMN_0004"],
                                               feature="CDS",
                                               logger=logger)[0]
        loc_index = coord_list[0]
        locus_tag = coord_list[4]
        strand = coord_list[2]
        coords = coord_list[1]
        seqid  = coord_list[5]
        self.assertEqual(loc_index, 0)
        self.assertEqual(locus_tag, "ECUMN_0004")
        self.assertEqual(strand, 1)
        self.assertEqual(coords, [3692, 4978])
        self.assertEqual(seqid, 'NC_011751.1')

    def test_pad_genbank_sequence(self):
        """Just tests the function, later we wll test that the same
        sequence is extracted padded and otherwise
        """
        padding_val = 50
        records = get_genbank_record(self.test_gb_file)
        old_seq = records[0].seq
        old_list = extract_coords_from_locus(record=records[0],
                                             locus_tag_list=["ECUMN_0004"],
                                             feature="CDS",
                                             logger=logger)
        old_start = old_list[0][1][0]
        coords, seq = pad_genbank_sequence(record=records[0],
                                           old_coords=old_list,
                                           padding=padding_val, logger=None)
        new_start = coords[0][1][0]
        self.assertEqual(old_start, new_start - padding_val)  # check index
        self.assertEqual(old_seq, seq[padding_val: -padding_val])  # full seq
        self.assertEqual(old_seq[:padding_val], seq[-padding_val:])  # 3' pad
        self.assertEqual(old_seq[-padding_val:], seq[:padding_val])  # 5' pad

    def test_strictly_increasing(self):
        self.assertTrue(strictly_increasing([1, 5, 5.5, 10, 10], dup_ok=True))
        with self.assertRaises(ValueError):
            strictly_increasing([1, 5, 5.5, 10, 10], dup_ok=False)
        self.assertFalse(strictly_increasing([1, 10, 5.5, 7, 9], dup_ok=False,
                                             verbose=False))

    def test_stitching(self):
        """  This is actually the thing needing the most testing, most likely
        """
        records = get_genbank_record(self.test_gb_file)
        clusters = parse_clustered_loci_file(self.test_loci_file,
                                             logger=logger)
        record = get_genbank_rec_from_multigb(recordID='NC_011751.1',
                                              genbank_records=records)
        coord_list = extract_coords_from_locus(record=record,
                                               locus_tag_list=clusters[0][1],
                                               feature="rRNA",
                                               logger=logger)
        stitched_record = stitch_together_target_regions(genome_sequence=\
                                                         record.seq,
                                                         coords=coord_list,
                                                         flanking="700:700",
                                                         within=50, minimum=50,
                                                         replace=False,
                                                         padding=0, # unused
                                                         logger=logger,
                                                         verbose=False)
        with open(self.test_cluster1, 'r') as ref:
            ref_rec = list(SeqIO.parse(ref, 'fasta'))[0]
        self.assertEqual(ref_rec.seq, stitched_record.seq)
        #TODO write test ccase with replacement

    def test_stitching_integration(self):
        """  Integration of several things
        """
        ex_padding = 1000  # an example padding amount
        records = get_genbank_record(self.test_gb_file)
        clusters = parse_clustered_loci_file(self.test_loci_file,
                                             logger=logger)
        record = get_genbank_rec_from_multigb(recordID='NC_011751.1',
                                              genbank_records=records)
        coord_list = extract_coords_from_locus(record=record,
                                               locus_tag_list=clusters[0][1],
                                               feature="rRNA",
                                               logger=logger)
        stitched_record = \
            stitch_together_target_regions(genome_sequence=\
                                           record.seq,
                                           coords=coord_list,
                                           flanking="700:700",
                                           within=50, minimum=50,
                                           replace=False,
                                           logger=logger,
                                           padding=ex_padding,
                                           verbose=False)
        with open(self.test_cluster1, 'r') as ref:
            ref_rec = list(SeqIO.parse(ref, 'fasta'))[0]
        self.assertEqual(ref_rec.seq, stitched_record.seq)
        padded_coords, padded_seq = pad_genbank_sequence(record=record,
                                                         old_coords=coord_list,
                                                         padding=ex_padding,
                                                         logger=None)

        # checks that the the sequence is properly padded
        self.assertEqual(record.seq,
                         padded_seq[ex_padding: -ex_padding])
        stitched_padded_record = \
            stitch_together_target_regions(genome_sequence=\
                                           padded_seq,
                                           coords=padded_coords,
                                           flanking="700:700",
                                           within=50, minimum=50,
                                           replace=False,
                                           logger=logger,
                                           verbose=False,
                                           padding=ex_padding,
                                           circular=True)
        # check the extracted sequences are still the same, which verifies the
        # coords were accuratly adjusted
        self.assertEqual(stitched_record.seq, stitched_padded_record.seq)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
