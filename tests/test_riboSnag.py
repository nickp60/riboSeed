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
import os
import unittest
import shutil

from Bio import SeqIO

from pyutilsnrw.utils3_5 import get_genbank_record, combine_contigs, md5

from riboSeed.riboSnag import parse_clustered_loci_file, \
    extract_coords_from_locus, strictly_increasing, \
    stitch_together_target_regions, get_genbank_rec_from_multigb,\
    pad_genbank_sequence, prepare_prank_cmd, prepare_mafft_cmd,\
    calc_Shannon_entropy, calc_entropy_msa, pure_python_plotting

from riboSeed.riboSnag import loci_cluster, locus


sys.dont_write_bytecode = True
logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class riboSnag_TestCase(unittest.TestCase):
    """ tests riboSnag
    """
    def setUp(self):
        self.curdir = os.getcwd()
        self.testdirname = os.path.join(os.path.dirname(__file__),
                                        "output_riboSnag_tests")
        self.test_snag_dir = os.path.join(os.path.dirname(__file__),
                                          str("references" + os.path.sep +
                                              "sample_snag_output"))
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
        self.test_prank_msa = os.path.join(os.path.dirname(__file__),
                                           str("references" + os.path.sep +
                                               "best_MSA.fasta.best.fas"))
        self.test_mafft_msa = os.path.join(os.path.dirname(__file__),
                                           str("references" + os.path.sep +
                                               "mafft_msa.fasta"))
        self.test_ent_seq = os.path.join(os.path.dirname(__file__),
                                         str("references" + os.path.sep +
                                             "test_ent_seq.txt"))
        self.plotting_script = os.path.join(os.path.dirname(__file__),
                                            str("references" + os.path.sep +
                                                "plotting_script.R"))
        self.samtools_exe = "samtools"
        self.prank_exe = "prank"
        self.mafft_exe = "mafft"
        self.maxDiff = None
        self.to_be_removed = []

    def test_parse_loci(self):
        """this checks the parsing of riboSelect ouput
        """
        clusters = parse_clustered_loci_file(self.test_loci_file,
                                             logger=logger)
        ref_seq_name = 'NC_011751.1'
        ref_loci_list = ['ECUMN_16S_6', 'ECUMN_23S_6', 'ECUMN_5S_7']
        test_locus_tags = [x.locus_tag for x in clusters[0].loci_list]
        self.assertEqual(clusters[0].sequence, ref_seq_name)
        self.assertEqual(test_locus_tags, ref_loci_list)

    # def test_genbank_match_id(self):
    #     """ tests get_genbank_rec_from_multigb
    #     """
    #     records = get_genbank_record(self.test_gb_file)
    #     record = get_genbank_rec_from_multigb(recordID='NC_011751.1',
    #                                           genbank_records=records)
    #     self.assertEqual(records[0].seq, record.seq)

    def test_extract_coords_from_locus(self):
        """todo: replace essentially unused get_genbank_record function
        """
        test_cluster = loci_cluster(index=0, sequence='NC_011751.1',
                                    loci_list=[locus(index=0,
                                                     sequence='NC_011751.1',
                                                     locus_tag="ECUMN_0004",
                                                     strand=None,
                                                     start_coord=None,
                                                     end_coord=None,
                                                     product=None)],
                                    padding=None,
                                    global_start_coord=None,
                                    global_end_coord=None,
                                    replace=False)
        records = get_genbank_record(self.test_gb_file)
        cluster_post_extract = extract_coords_from_locus(cluster=test_cluster,
                                                         record=records[0],
                                                         feature="CDS",
                                                         logger=logger)
        self.assertEqual(cluster_post_extract.loci_list[0].index, 0)
        self.assertEqual(
            cluster_post_extract.loci_list[0].locus_tag, "ECUMN_0004")
        self.assertEqual(
            cluster_post_extract.loci_list[0].strand, 1)
        self.assertEqual(
            cluster_post_extract.loci_list[0].start_coord, 3692)
        self.assertEqual(
            cluster_post_extract.loci_list[0].end_coord, 4978)
        self.assertEqual(
            cluster_post_extract.loci_list[0].sequence, 'NC_011751.1')

    # def test_pad_genbank_sequence(self):
    #     """Just tests the function, later we wll test that the same
    #     sequence is extracted padded and otherwise
    #     """
    #     padding_val = 50
    #     records = get_genbank_record(self.test_gb_file)
    #     old_seq = records[0].seq
    #     old_list = extract_coords_from_locus(record=records[0],
    #                                          locus_tag_list=["ECUMN_0004"],
    #                                          feature="CDS",
    #                                          logger=logger)
    #     old_start = old_list[0][1][0]
    #     coords, seq = pad_genbank_sequence(record=records[0],
    #                                        old_coords=old_list,
    #                                        padding=padding_val, logger=None)
    #     new_start = coords[0][1][0]
    #     self.assertEqual(old_start, new_start - padding_val)  # check index
    #     self.assertEqual(old_seq, seq[padding_val: -padding_val])  # full seq
    #     self.assertEqual(old_seq[:padding_val], seq[-padding_val:])  # 3' pad
    #     self.assertEqual(old_seq[-padding_val:], seq[:padding_val])  # 5' pad

    # def test_strictly_increasing(self):
    #     """ I pulled this largely from SO, and it should check to
    #     see if items are in an increasing order, with and without duplicates
    #     """
    #     self.assertTrue(strictly_increasing([1, 5, 5.5, 10, 10], dup_ok=True))
    #     with self.assertRaises(ValueError):
    #         strictly_increasing([1, 5, 5.5, 10, 10], dup_ok=False)
    #     self.assertFalse(strictly_increasing([1, 10, 5.5, 7, 9], dup_ok=False,
    #                                          verbose=False))

    # def test_stitching(self):
    #     """  This is actually the thing needing the most testing, most likely
    #     """
    #     records = get_genbank_record(self.test_gb_file)
    #     clusters = parse_clustered_loci_file(self.test_loci_file,
    #                                          logger=logger)
    #     record = get_genbank_rec_from_multigb(recordID='NC_011751.1',
    #                                           genbank_records=records)
    #     coord_list = extract_coords_from_locus(record=record,
    #                                            locus_tag_list=clusters[0][1],
    #                                            feature="rRNA",
    #                                            logger=logger)
    #     stitched_record = \
    #         stitch_together_target_regions(genome_sequence=record.seq,
    #                                        coords=coord_list,
    #                                        flanking="700:700",
    #                                        within=50, minimum=50,
    #                                        replace=False,
    #                                        padding=0,  # unused
    #                                        logger=logger,
    #                                        verbose=False)
    #     with open(self.test_cluster1, 'r') as ref:
    #         ref_rec = list(SeqIO.parse(ref, 'fasta'))[0]
    #     self.assertEqual(ref_rec.seq, stitched_record.seq)
    #     # If reimplementing replacement write test cases here

    # def test_stitching_integration(self):
    #     """  Integration of several things
    #     """
    #     ex_padding = 1000  # an example padding amount
    #     records = get_genbank_record(self.test_gb_file)
    #     clusters = parse_clustered_loci_file(self.test_loci_file,
    #                                          logger=logger)
    #     record = get_genbank_rec_from_multigb(recordID='NC_011751.1',
    #                                           genbank_records=records)
    #     coord_list = extract_coords_from_locus(record=record,
    #                                            locus_tag_list=clusters[0][1],
    #                                            feature="rRNA",
    #                                            logger=logger)
    #     stitched_record = \
    #         stitch_together_target_regions(genome_sequence=record.seq,
    #                                        coords=coord_list,
    #                                        flanking="700:700",
    #                                        within=50, minimum=50,
    #                                        replace=False,
    #                                        logger=logger,
    #                                        padding=ex_padding,
    #                                        verbose=False)
    #     with open(self.test_cluster1, 'r') as ref:
    #         ref_rec = list(SeqIO.parse(ref, 'fasta'))[0]
    #     self.assertEqual(ref_rec.seq, stitched_record.seq)
    #     padded_coords, padded_seq = pad_genbank_sequence(record=record,
    #                                                      old_coords=coord_list,
    #                                                      padding=ex_padding,
    #                                                      logger=None)

    #     # checks that the the sequence is properly padded
    #     self.assertEqual(record.seq,
    #                      padded_seq[ex_padding: -ex_padding])
    #     stitched_padded_record = \
    #         stitch_together_target_regions(genome_sequence=padded_seq,
    #                                        coords=padded_coords,
    #                                        flanking="700:700",
    #                                        within=50, minimum=50,
    #                                        replace=False,
    #                                        logger=logger,
    #                                        verbose=False,
    #                                        padding=ex_padding,
    #                                        circular=True)
    #     # check the extracted sequences are still the same, which verifies the
    #     # coords were accuratly adjusted
    #     self.assertEqual(stitched_record.seq, stitched_padded_record.seq)

    # def test_calc_Shannon(self):
    #     """ test calc_shannon_entropy;
    #     may need to recalculate
    #     """
    #     test_matrix = [
    #         ["A", "A", "A", "A", "A"],
    #         ["A", "A", "A", "A", "G"],
    #         ["A", "A", "A", "G", "G"],
    #         ["A", "A", "G", "G", "G"],
    #         ["A", "T", "G", "G", "G"],
    #         ["T", "T", "G", "G", "G"],
    #         ["A", "C", "G", "T", "-"]
    #     ]
    #     entropies = calc_Shannon_entropy(test_matrix)
    #     print(entropies)
    #     theoretical = [-0.0, 0.500402, 0.673012, 0.673012,
    #                    0.950271, 0.673012, 1.609438]
    #     for index, value in enumerate(entropies):
    #         self.assertEqual(round(value, 6), theoretical[index])

    # def test_calc_entropy_from_msa(self):
    #     """
    #     """
    #     # calc_Shannon_entropy(msa_path=self.test_prank_msa)
    #     values = calc_entropy_msa(msa_path=self.test_mafft_msa)
    #     with open(self.test_ent_seq, "r") as efile:
    #         for index, line in enumerate(efile):
    #             print(line)
    #             print(line.strip())
    #             self.assertEqual(round(float(line.strip()), 7),
    #                              round(values[index], 7))

    # def test_prepare_msa_cmds(self):
    #     """ test cases of prank and mafft command creation
    #     """
    #     for dir in [self.testdirname, self.test_snag_dir]:
    #         try:
    #             os.makedirs(self.testdirname)
    #         except:
    #             print("using existing {0} directory".format(dir))
    #             pass
    #     unaligned_seqs = combine_contigs(contigs_dir=self.test_snag_dir,
    #                                      pattern="*",
    #                                      contigs_name="riboSnag_unaligned",
    #                                      ext=".fasta", verbose=False,
    #                                      logger=logger)
    #     prank_cmd_1, results_path = prepare_prank_cmd(
    #         outdir=self.testdirname,
    #         outfile_name="best_MSA.fasta",
    #         combined_fastas=unaligned_seqs,
    #         prank_exe=self.prank_exe,
    #         add_args="-t=sometree",
    #         clobber=False, logger=None)
    #     mafft_cmd_1, results_path = prepare_mafft_cmd(
    #         outdir=self.testdirname,
    #         outfile_name="best_MSA.fasta",
    #         combined_fastas=unaligned_seqs,
    #         mafft_exe=self.mafft_exe,
    #         add_args="-t=sometree",
    #         clobber=False, logger=None)
    #     idealprank = "prank -t=sometree -d={0} -o={1}".format(
    #         unaligned_seqs, os.path.join(self.testdirname,
    #                                      "best_MSA.fasta"))
    #     idealmafft = "mafft -t=sometree {0} > {1}".format(
    #         unaligned_seqs, os.path.join(self.testdirname,
    #                                      "best_MSA.fasta"))
    #     self.assertEqual(idealprank, prank_cmd_1)
    #     self.assertEqual(idealmafft, mafft_cmd_1)
    #     self.to_be_removed.append(unaligned_seqs)

    # def test_pure_python_plotting(self):
    #     """
    #     """
    #     entropies = []
    #     with open(self.test_ent_seq, 'r') as entfile:
    #         for line in entfile:
    #             entropies.append(float(line.strip()))
    #     pure_python_plotting(data=entropies, script_name="testplot.R",
    #                          outdir=self.testdirname, name="Sample",
    #                          outfile_prefix="test_generated_entropy_plot",
    #                          DEBUG=True)

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            os.unlink(filename)

if __name__ == '__main__':
    unittest.main()
