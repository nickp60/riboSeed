# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas
"""
import copy
import sys
import logging
import os
import unittest
import shutil
import time



from Bio import SeqIO


from riboSeed.shared_methods import parse_clustered_loci_file, \
    pad_genbank_sequence, get_genbank_record, combine_contigs, md5, \
    extract_coords_from_locus, add_gb_seqrecords_to_cluster_list

from riboSeed.riboSnag import stitch_together_target_regions, \
    prepare_prank_cmd, prepare_mafft_cmd,\
    calc_Shannon_entropy, calc_entropy_msa,\
    annotate_msa_conensus, plot_scatter_with_anno, get_all_kmers,\
    profile_kmer_occurances, plot_pairwise_least_squares, make_msa, \
    check_loci_file_not_genbank

from riboSeed.classes import LociCluster, Locus


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
        self.test_cluster2 = os.path.join(os.path.dirname(__file__),
                                          str("references" + os.path.sep +
                                              'cluster2.fasta'))
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
        self.maxDiff = 1000
        self.startTime = time.time()
        self.to_be_removed = []
        if not os.path.exists(self.testdirname):
            os.makedirs(self.testdirname, exist_ok=True)

    def test_check_loci_not_gb(self):
        """ check that the loci file is not a genbank file
        """
        with self.assertRaises(FileNotFoundError):
            check_loci_file_not_genbank(self.test_gb_file)

    def test_parse_loci(self):
        """this checks the parsing of riboSelect ouput
        """
        clusters = parse_clustered_loci_file(
            filepath=self.test_loci_file,
            gb_filepath=self.test_gb_file,
            output_root=self.testdirname,
            padding=100,
            circular=True,
            logger=logger)
        ref_seq_name = 'NC_011751.1'
        ref_loci_list = ['ECUMN_16S_6', 'ECUMN_23S_6', 'ECUMN_5S_7']
        test_locus_tags = [x.locus_tag for x in clusters[0].loci_list]
        self.assertEqual(clusters[0].sequence_id, ref_seq_name)
        self.assertEqual(test_locus_tags, ref_loci_list)

    # def test_genbank_match_id(self):
    #     """ tests get_genbank_rec_from_multigb
    #     """
    #     records = get_genbank_record(self.test_gb_file)
    #     record = get_genbank_rec_from_multigb(recordID='NC_011751.1',
    #                                           genbank_records=records)
    #     self.assertEqual(records[0].seq, record.seq)
    #     with self.assertRaises(ValueError):
    #         get_genbank_rec_from_multigb(
    #             recordID='NC_011751.X', genbank_records=records)

    def test_extract_coords_from_locus(self):
        """todo: replace essentially unused get_genbank_record function
        """
        test_cluster = LociCluster(
            sequence_id='NC_011751.1',
            loci_list=[Locus(index=0,
                             sequence_id='NC_011751.1',
                             locus_tag="ECUMN_0004",
                             strand=None,
                             start_coord=None,
                             end_coord=None,
                             product=None)],
            padding=None,
            global_start_coord=None,
            global_end_coord=None)
        test_cluster.seq_record = get_genbank_record(self.test_gb_file)[0]
        extract_coords_from_locus(
            cluster=test_cluster,
            feature="CDS",
            logger=logger)
        cluster_post_extract = test_cluster
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
            cluster_post_extract.loci_list[0].sequence_id, 'NC_011751.1')

    def test_pad_genbank_sequence(self):
        """Just tests the function, later we wll test that the same
        sequence is extracted padded and otherwise
        """
        padding_val = 50
        test_cluster = LociCluster(
            sequence_id='NC_011751.1',
            loci_list=[Locus(index=0,
                             sequence_id='NC_011751.1',
                             locus_tag="ECUMN_0004",
                             strand=None,
                             start_coord=None,
                             end_coord=None,
                             product=None)],
            padding=padding_val,
            global_start_coord=None,
            global_end_coord=None)
        test_cluster.seq_record = get_genbank_record(self.test_gb_file)[0]
        extract_coords_from_locus(cluster=test_cluster,
                                  feature="CDS",
                                  logger=logger)
        old_seq = test_cluster.seq_record
        # test excessive padding
        cluster_too_much_padding = copy.deepcopy(test_cluster)
        cluster_too_much_padding.padding = 10000000
        # with self.assertRaises(ValueError):
        #     pad_genbank_sequence(cluster_too_much_padding,
        #                          logger=logger)

        # test sucessful execution
        padded_cluster = pad_genbank_sequence(test_cluster,
                                              logger=logger)
        self.assertEqual(str(old_seq.seq),
                         str(padded_cluster.seq_record.seq[padding_val:
                                                          - padding_val]))

    def test_stitching(self):
        """  This is actually the thing needing the most testing, most likely
        """
        padding_val = 1000
        clusters = parse_clustered_loci_file(
            filepath=self.test_loci_file,
            gb_filepath=self.test_gb_file,
            output_root=self.testdirname,
            padding=padding_val,
            circular=True,
            logger=logger)
        clusters = add_gb_seqrecords_to_cluster_list(
            cluster_list=clusters,
            gb_filepath=self.test_gb_file)
        cluster = clusters[0]  # this should be reverse complimented
        extract_coords_from_locus(
            cluster=cluster,
            feature="rRNA",
            logger=logger)
        stitched_cluster = stitch_together_target_regions(
            cluster=cluster,
            flanking=1000,
            logger=logger,
            circular=False,
            revcomp=True)
        stitched_cluster_circular = stitch_together_target_regions(
            cluster=cluster,
            flanking=1000,
            logger=logger,
            circular=True,
            revcomp=True)
        # # fail with invalid flanking arg
        # with self.assertRaises(ValueError):
        #     stitch_together_target_regions(
        #         cluster=cluster,
        #         flanking="1000:1000",
        #         logger=logger,
        #         circular=False,
        #         revcomp=False)
        # # fail with exceeded minimum
        # with self.assertRaises(ValueError):
        #     stitch_together_target_regions(
        #         cluster=cluster,
        #         flanking="1000",
        #         within=5000, minimum=1500,
        #         replace=True,
        #         logger=logger,
        #         circular=False,
        #         revcomp=False)
        with open(self.test_cluster2, 'r') as ref:
            ref_rec = list(SeqIO.parse(ref, 'fasta'))[0]
        self.assertEqual(str(ref_rec.seq),
                         str(stitched_cluster.extractedSeqRecord.seq))
        # ensure that coords are named identically whether genome is
        #  treated as circular
        self.assertEqual(stitched_cluster.extractedSeqRecord.id,
                         stitched_cluster_circular.extractedSeqRecord.id)

    def test_calc_Shannon(self):
        """ test calc_shannon_entropy;
        may need to recalculate
        """
        test_matrix = [
            ["A", "A", "A", "A", "A"],
            ["A", "A", "A", "A", "G"],
            ["A", "A", "A", "G", "G"],
            ["A", "A", "G", "G", "G"],
            ["A", "T", "G", "G", "G"],
            ["T", "T", "G", "G", "G"],
            ["A", "C", "G", "T", "-"]
        ]
        entropies = calc_Shannon_entropy(matrix=test_matrix)
        theoretical = [-0.0, 0.500402, 0.673012, 0.673012,
                       0.950271, 0.673012, 1.609438]
        for index, value in enumerate(entropies):
            self.assertEqual(round(value, 6), theoretical[index])

    def test_calc_entropy_from_msa(self):
        """
        """
        # calc_Shannon_entropy(msa_path=self.test_prank_msa)
        seq_entropy, names, tseq = calc_entropy_msa(
            msa_path=self.test_mafft_msa)
        with open(self.test_ent_seq, "r") as efile:
            for index, line in enumerate(efile):
                self.assertEqual(round(float(line.strip()), 7),
                                 round(seq_entropy[index], 7))

    def test_prepare_msa_cmds(self):
        """ test cases of prank and mafft command creation
        """
        for dr in [self.test_snag_dir]:
            try:
                os.makedirs(dr, exist_ok=True)
            except:
                print("\nusing existing {0} directory".format(dir))
                pass
        unaligned_seqs = combine_contigs(contigs_dir=self.test_snag_dir,
                                         pattern="*",
                                         contigs_name="riboSnag_unaligned",
                                         ext=".fasta", verbose=False,
                                         logger=logger)
        prank_cmd_1, results_path = prepare_prank_cmd(
            outdir=self.testdirname,
            outfile_name="best_MSA",
            combined_fastas=unaligned_seqs,
            prank_exe=self.prank_exe,
            add_args="-t=sometree",
            logger=logger)
        mafft_cmd_1, results_path = prepare_mafft_cmd(
            outdir=self.testdirname,
            outfile_name="best_MSA",
            combined_fastas=unaligned_seqs,
            mafft_exe=self.mafft_exe,
            add_args="-t=sometree",
            logger=logger)
        idealprank = "prank -t=sometree -d={0} -o={1}".format(
            unaligned_seqs, os.path.join(self.testdirname,
                                         "best_MSA"))
        idealmafft = "mafft -t=sometree {0} > {1}".format(
            unaligned_seqs, os.path.join(self.testdirname,
                                         "best_MSA.best.fas"))
        self.assertEqual(idealprank, prank_cmd_1)
        self.assertEqual(idealmafft, mafft_cmd_1)
        self.to_be_removed.append(unaligned_seqs)

    @unittest.skipIf(shutil.which("mafft") is None or
                     shutil.which("prank") is None,
                     "maft or prank  executables not found, skipping."+
                     "If this isnt an error from travis deployment, you probably "+
                     "should install it")
    def test_make_msa(self):
        """
        """
        unaligned_seqs = combine_contigs(contigs_dir=self.test_snag_dir,
                                         pattern="*",
                                         contigs_name="riboSnag_unaligned",
                                         ext=".fasta", verbose=False,
                                         logger=logger)
        # with mafft
        msa_cmd1, results_path1 = make_msa(msa_tool='mafft',
                                           unaligned_seqs=unaligned_seqs,
                                           args="-t=sometree",
                                           prank_exe=self.prank_exe,
                                           mafft_exe=self.mafft_exe,
                                           outdir=self.testdirname,
                                           logger=logger)
        # with prank
        msa_cmd2, results_path2 = make_msa(msa_tool='prank',
                                           unaligned_seqs=unaligned_seqs,
                                           args="-t=sometree",
                                           prank_exe=self.prank_exe,
                                           mafft_exe=self.mafft_exe,
                                           outdir=self.testdirname,
                                           logger=logger)
        idealprank = "prank -t=sometree -d={0} -o={1}".format(
            unaligned_seqs, os.path.join(self.testdirname,
                                         "best_MSA"))
        idealmafft = "mafft -t=sometree {0} > {1}".format(
            unaligned_seqs, os.path.join(self.testdirname,
                                         "best_MSA.best.fas"))
        self.assertEqual(results_path2, results_path1)
        self.assertEqual(idealprank, msa_cmd2)
        self.assertEqual(idealmafft, msa_cmd1)
        self.to_be_removed.append(unaligned_seqs)

    @unittest.skipIf(shutil.which("barrnap") is None,
                     "barnnap executable not found, skipping."+
                     "If this isnt an error from travis deployment, you probably "+
                     "should install it")
    def test_msa_consensus(self):
        """ calculate entropy, annotate consenesus,
        and ensure barrnap is wired up properly
        """
        seq_entropy, names, tseq = calc_entropy_msa(
            msa_path=self.test_mafft_msa)
        with open(self.test_ent_seq, "r") as efile:
            for index, line in enumerate(efile):
                self.assertEqual(round(float(line.strip()), 7),
                                 round(seq_entropy[index], 7))
        gff, conseq, named_coords = annotate_msa_conensus(
            tseq_array=tseq,
            seq_file=os.path.join(self.testdirname,
                                  "consensus_sample.fasta"),
            pattern='product=(.+?)$',
            barrnap_exe="barrnap",
            kingdom="bact",
            collapseNs=False,
            countdashcov=True,
            excludedash=False,
            logger=logger)
        # check first entry has 9 fields
        self.assertEqual(len(gff[1]), 9)
        # check start coord is still 362
        self.assertEqual(gff[1][4], '362')
        # check name
        self.assertEqual(named_coords[0][0], '5S ribosomal RNA')

    @unittest.skipIf(shutil.which("barrnap") is None,
                     "barnnap executable not found, skipping."+
                     "If this isnt an error from travis deployment, you probably "+
                     "should install it")
    def test_scatter_plotting(self):
        """
        """
        seq_entropy, names, tseq = calc_entropy_msa(
            msa_path=self.test_mafft_msa)
        with open(self.test_ent_seq, "r") as efile:
            for index, line in enumerate(efile):
                self.assertEqual(round(float(line.strip()), 7),
                                 round(seq_entropy[index], 7))
        gff, conseq, named_coords = annotate_msa_conensus(
            tseq_array=tseq,
            seq_file=os.path.join(self.testdirname,
                                  "consensus.fasta"),
            pattern='product=(.+?)$',
            barrnap_exe="barrnap",
            kingdom="bact",
            collapseNs=False,
            countdashcov=True,
            excludedash=False,
            logger=logger)
        plot_scatter_with_anno(data=seq_entropy,
                               names=["Position", "Entropy"],
                               title="Shannon Entropy by Position",
                               consensus_cov=conseq,
                               anno_list=named_coords,
                               output_prefix=os.path.join(self.testdirname,
                                                          "entropy_plot"))
        #  and without dashes
        gff_nodash, conseq_nodash, named_coords_nodash = annotate_msa_conensus(
            tseq_array=tseq,
            seq_file=os.path.join(self.testdirname,
                                  "consensus_nodash.fasta"),
            pattern='product=(.+?)$',
            barrnap_exe="barrnap",
            kingdom="bact",
            collapseNs=True,
            countdashcov=False,
            excludedash=False,
            logger=logger)
        plot_scatter_with_anno(
            data=seq_entropy,
            names=["Position", "Entropy"],
            title="Shannon Entropy by Position No N's",
            consensus_cov=conseq_nodash,
            anno_list=named_coords_nodash,
            output_prefix=os.path.join(self.testdirname, "entropy_plot_noN"))

    def test_get_all_kmers(self):
        """
        """
        string="ACGTCGACGACAGCTGAGCTGTCGTCTGCTGCGCTTA-T-TGCGACGTTACG"
        all_kmers = sorted(get_all_kmers(alph="ATCG-", length=2))
        ref_all_mers = ['--', '-A', '-C', '-G', '-T', 'A-', 'AA', 'AC', 'AG',
                        'AT', 'C-', 'CA', 'CC', 'CG', 'CT', 'G-', 'GA', 'GC',
                        'GG', 'GT', 'T-', 'TA', 'TC', 'TG', 'TT']
        self.assertEqual(all_kmers, ref_all_mers)

    def test_profile_kmer_occurances(self):
        """
        """
        with open(self.test_mafft_msa, 'r') as resfile:
            kmer_seqs = list(SeqIO.parse(resfile, "fasta"))
        occurances, seqnames = profile_kmer_occurances(
            rec_list=kmer_seqs,
            k=4,
            logger=logger)
        self.assertEqual(sorted(['NZ_CP017149.1_726410..727025',
                                 'NZ_CP017149.1_4787826..4788441_RC',
                                 'NZ_CP017149.1_5275421..5276036_RC',
                                 'NZ_CP017149.1_6050886..6051501_RC']),
                         sorted(seqnames))
        # TODO make test for the occurances

    def test_plot_least_squares(self):
        """the whole sls_matrix
        """
        with open(self.test_mafft_msa, 'r') as resfile:
            kmer_seqs = list(SeqIO.parse(resfile, "fasta"))
        occurances, seqnames = profile_kmer_occurances(
            rec_list=kmer_seqs,
            k=4,
            logger=logger)
        df = plot_pairwise_least_squares(
            counts=occurances,
            names_list=seqnames,
            output_prefix=os.path.join(
                self.testdirname, "sum_least_squares"))
        sls_matrix = [['NZ_CP017149.1_726410..727025',
                       'NZ_CP017149.1_726410..727025', 0],
                      ['NZ_CP017149.1_726410..727025',
                       'NZ_CP017149.1_4787826..4788441_RC', 512],
                      ['NZ_CP017149.1_726410..727025',
                       'NZ_CP017149.1_5275421..5276036_RC', 448],
                      ['NZ_CP017149.1_726410..727025',
                       'NZ_CP017149.1_6050886..6051501_RC', 450],
                      ['NZ_CP017149.1_4787826..4788441_RC',
                       'NZ_CP017149.1_726410..727025', 512],
                      ['NZ_CP017149.1_4787826..4788441_RC',
                       'NZ_CP017149.1_4787826..4788441_RC', 0],
                      ['NZ_CP017149.1_4787826..4788441_RC',
                       'NZ_CP017149.1_5275421..5276036_RC', 362],
                      ['NZ_CP017149.1_4787826..4788441_RC',
                       'NZ_CP017149.1_6050886..6051501_RC', 380],
                      ['NZ_CP017149.1_5275421..5276036_RC',
                       'NZ_CP017149.1_726410..727025', 448],
                      ['NZ_CP017149.1_5275421..5276036_RC',
                       'NZ_CP017149.1_4787826..4788441_RC', 362],
                      ['NZ_CP017149.1_5275421..5276036_RC',
                       'NZ_CP017149.1_5275421..5276036_RC', 0],
                      ['NZ_CP017149.1_5275421..5276036_RC',
                       'NZ_CP017149.1_6050886..6051501_RC', 392],
                      ['NZ_CP017149.1_6050886..6051501_RC',
                       'NZ_CP017149.1_726410..727025', 450],
                      ['NZ_CP017149.1_6050886..6051501_RC',
                       'NZ_CP017149.1_4787826..4788441_RC', 380],
                      ['NZ_CP017149.1_6050886..6051501_RC',
                       'NZ_CP017149.1_5275421..5276036_RC', 392],
                      ['NZ_CP017149.1_6050886..6051501_RC',
                       'NZ_CP017149.1_6050886..6051501_RC', 0]]
        for index in range(0, len(sls_matrix)):
            self.assertEqual(df.as_matrix().tolist()[index],
                             sls_matrix[index])

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            os.unlink(filename)
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))

if __name__ == '__main__':
    unittest.main()
