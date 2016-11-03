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
    extract_coords_from_locus, \
    stitch_together_target_regions, get_genbank_rec_from_multigb,\
    pad_genbank_sequence, prepare_prank_cmd, prepare_mafft_cmd,\
    calc_Shannon_entropy, calc_entropy_msa,\
    annotate_msa_conensus, plot_scatter_with_anno, get_all_kmers,\
    profile_kmer_occurances, plot_pairwise_least_squares

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
        self.to_be_removed = []

    def test_parse_loci(self):
        """this checks the parsing of riboSelect ouput
        """
        # error parsing
        with self.assertRaises(ValueError):
            parse_clustered_loci_file(
                filepath=str(self.test_loci_file + "png"),
                gb_filepath=self.test_gb_file,
                padding=100,
                circular=True,
                logger=logger)
        # error wrong files provided
        with self.assertRaises(FileNotFoundError):
            parse_clustered_loci_file(
                gb_filepath=self.test_loci_file,
                filepath=self.test_gb_file,
                padding=100,
                circular=True,
                logger=logger)
        clusters = parse_clustered_loci_file(
            filepath=self.test_loci_file,
            gb_filepath=self.test_gb_file,
            padding=100,
            circular=True,
            logger=logger)
        ref_seq_name = 'NC_011751.1'
        ref_loci_list = ['ECUMN_16S_6', 'ECUMN_23S_6', 'ECUMN_5S_7']
        test_locus_tags = [x.locus_tag for x in clusters[0].loci_list]
        self.assertEqual(clusters[0].sequence, ref_seq_name)
        self.assertEqual(test_locus_tags, ref_loci_list)

    def test_genbank_match_id(self):
        """ tests get_genbank_rec_from_multigb
        """
        records = get_genbank_record(self.test_gb_file)
        record = get_genbank_rec_from_multigb(recordID='NC_011751.1',
                                              genbank_records=records)
        self.assertEqual(records[0].seq, record.seq)
        with self.assertRaises(ValueError):
            get_genbank_rec_from_multigb(
                recordID='NC_011751.X', genbank_records=records)
    def test_extract_coords_from_locus(self):
        """todo: replace essentially unused get_genbank_record function
        """
        test_cluster = loci_cluster(
            index=0, sequence='NC_011751.1',
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
        test_cluster.SeqRecord = get_genbank_record(self.test_gb_file)[0]
        cluster_post_extract = extract_coords_from_locus(
            cluster=test_cluster,
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

    def test_pad_genbank_sequence(self):
        """Just tests the function, later we wll test that the same
        sequence is extracted padded and otherwise
        """
        padding_val = 50
        test_cluster = loci_cluster(
            index=0, sequence='NC_011751.1',
            loci_list=[locus(index=0,
                             sequence='NC_011751.1',
                             locus_tag="ECUMN_0004",
                             strand=None,
                             start_coord=None,
                             end_coord=None,
                             product=None)],
            padding=padding_val,
            global_start_coord=None,
            global_end_coord=None,
            replace=False)
        test_cluster.SeqRecord = get_genbank_record(self.test_gb_file)[0]
        cluster_post_extract = extract_coords_from_locus(cluster=test_cluster,
                                                         feature="CDS",
                                                         logger=logger)
        old_seq = cluster_post_extract.SeqRecord
        # old_list = [[x.start_coord, x.end_coord] for
        #             x in cluster_post_extract.loci_list]
        padded_cluster = pad_genbank_sequence(cluster_post_extract,
                                              logger=logger)
        self.assertEqual(str(old_seq.seq),
                         str(padded_cluster.SeqRecord.seq[padding_val:
                                                          -padding_val]))


    def test_stitching(self):
        """  This is actually the thing needing the most testing, most likely
        """
        padding_val = 1000
        clusters = parse_clustered_loci_file(
            filepath=self.test_loci_file,
            gb_filepath=self.test_gb_file,
            padding=padding_val,
            circular=True,
            logger=logger)
        cluster = clusters[0]
        cluster_post_extract = extract_coords_from_locus(
            cluster=cluster,
            feature="rRNA",
            logger=logger)
        stitched_cluster = stitch_together_target_regions(
            cluster=cluster_post_extract,
            flanking="1000",
            within=50, minimum=50,
            replace=False,
            logger=logger,
            verbose=False,
            circular=False,
            revcomp=False)
        with open(self.test_cluster2, 'r') as ref:
            ref_rec = list(SeqIO.parse(ref, 'fasta'))[0]
        self.assertEqual(str(ref_rec.seq),
                         str(stitched_cluster.extractedSeqRecord.seq))
        # If reimplementing replacement write test cases here

    # def test_stitching_integration(self):
    #     """  Integration of several things
    #     """
    #     ex_padding = 1000  # an example padding amount
    #     clusters = parse_clustered_loci_file(
    #         filepath=self.test_loci_file,
    #         gb_filepath=self.test_gb_file,
    #         padding=ex_padding,
    #         circular=False,
    #         logger=logger)
    #     cluster = clusters[0]
    #     cluster_post_extract = extract_coords_from_locus(
    #         cluster=cluster,
    #         feature="rRNA",
    #         logger=logger)
    #     stitched_cluster = stitch_together_target_regions(
    #         cluster=cluster_post_extract,
    #         flanking="1000:1000",
    #         within=50, minimum=50,
    #         replace=False,
    #         logger=logger,
    #         verbose=False,
    #         circular=False,
    #         revcomp=False)
    #     with open(self.test_cluster2, 'r') as ref:
    #         ref_rec = list(SeqIO.parse(ref, 'fasta'))[0]
    #     self.assertEqual(str(ref_rec.seq),
    #                      str(stitched_cluster.extractedSeqRecord.seq))

    #     padded_cluster = pad_genbank_sequence(stiched_cluster,
    #                                           logger=logger)

    #     # checks that the the sequence is properly padded
    #     records = get_genbank_record(self.test_gb_file)
    #     record = get_genbank_rec_from_multigb(recordID='NC_011751.1',
    #                                           genbank_records=records)

    #     self.assertEqual(record.seq,
    #                      padded_cluster.SeqRecord.seq[ex_padding: -ex_padding])
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
        entropies = calc_Shannon_entropy(test_matrix)
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
        for dir in [self.testdirname, self.test_snag_dir]:
            try:
                os.makedirs(self.testdirname)
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
            outfile_name="best_MSA.fasta",
            combined_fastas=unaligned_seqs,
            prank_exe=self.prank_exe,
            add_args="-t=sometree",
            clobber=False, logger=None)
        mafft_cmd_1, results_path = prepare_mafft_cmd(
            outdir=self.testdirname,
            outfile_name="best_MSA.fasta",
            combined_fastas=unaligned_seqs,
            mafft_exe=self.mafft_exe,
            add_args="-t=sometree",
            clobber=False, logger=None)
        idealprank = "prank -t=sometree -d={0} -o={1}".format(
            unaligned_seqs, os.path.join(self.testdirname,
                                         "best_MSA.fasta"))
        idealmafft = "mafft -t=sometree {0} > {1}".format(
            unaligned_seqs, os.path.join(self.testdirname,
                                         "best_MSA.fasta"))
        self.assertEqual(idealprank, prank_cmd_1)
        self.assertEqual(idealmafft, mafft_cmd_1)
        self.to_be_removed.append(unaligned_seqs)

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            os.unlink(filename)

if __name__ == '__main__':
    unittest.main()


# this is kept around for sentimental reasons
# def pure_python_plotting(data, outdir, script_name="plot.R", name="",
#                          outfile_prefix="entropy_plot", DEBUG=True):
#     """
#     """
#     if not os.path.exists(outdir):
#         raise ValueError("Output directory not found")
#     datafile = os.path.join(outdir, str(outfile_prefix + ".csv"))
#     imgfile = os.path.join(outdir, str(outfile_prefix + ".pdf"))
#     with open(datafile, "w") as df:
#         df.write(str("position, entropy\n"))
#         for en, i in enumerate(data):
#             df.write("{0}, {1}\n".format(str(en), str(i)))
#     rcmds = ["# Generated by riboSnag.py on {0}".format(time.asctime()),
#              str("data <- read.csv('{0}', stringsAsFactors = F)").format(
#                  datafile),
#              str("pdf(file='{0}', width=6, height=4)").format(imgfile),
#              "plot(data, cex=.75, pch=18, axes=FALSE, xaxs='i', yaxs='i'," +
#              " ylim=c(-0.1, 2), xlab='', ylab='')",
#              str("mtext('{0}', side=3, line=0.5, cex.lab=1,las=1, " +
#                  "col='#34282C')").format(name),
#              "axis(2, col='darkgrey')",
#              "axis(1, col='darkgrey')",
#              "title('Shannon Entropy', xlab='Position (bp)', ylab='Entropy')",
#              "dev.off()"
#              ]
#     with open(os.path.join(os.getcwd(), script_name), "w") as rf:
#         for cmd in rcmds:
#             rf.write(cmd + "\n")
#     subprocess.run("Rscript {0}".format(script_name),
#                    shell=sys.platform != 'win32',
#                    check=True)
#     if not DEBUG:
#         os.remove(os.path.join(os.getcwd(), datafile))
#         os.remove(os.path.join(os.getcwd(), script_name))


#     def test_pure_python_plotting(self):
#         """I kept this around for historical reasons; function is below
#         """
#         entropies = []
#         with open(self.test_ent_seq, 'r') as entfile:
#             for line in entfile:
#                 entropies.append(float(line.strip()))
#         pure_python_plotting(data=entropies, script_name="testplot.R",
#                              outdir=self.testdirname, name="Sample",
#                              outfile_prefix="test_generated_entropy_plot",
#                              DEBUG=True)
