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
import multiprocessing
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import Namespace

# I hate this line but it works :(
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "riboSeed"))


from pyutilsnrw.utils3_5 import md5, file_len, copy_file


from riboSeed.riboSeed2 import SeedGenome, ngsLib,  LociMapping, \
    map_to_genome_ref_smalt, add_coords_to_clusters, partition_mapping, \
     convert_bams_to_fastq_cmds, \
    generate_spades_cmd, estimate_distances_smalt, run_final_assemblies


from riboSeed.riboSnag import parse_clustered_loci_file, \
    extract_coords_from_locus, \
    stitch_together_target_regions, get_genbank_rec_from_multigb, \
    pad_genbank_sequence, prepare_prank_cmd, prepare_mafft_cmd, \
    calc_Shannon_entropy, calc_entropy_msa, \
    annotate_msa_conensus, plot_scatter_with_anno, \
    profile_kmer_occurances, plot_pairwise_least_squares, make_msa


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
        self.ref_fasta = os.path.join(self.test_dir,
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
        self.test_loci_file = os.path.join(os.path.dirname(__file__),
                                           str("references" + os.path.sep +
                                               'grouped_loci_reference.txt'))
        self.args = Namespace(skip_contol=False, kmers="21,33,55,77,99",
                              spades_exe="spades.py", quast_exe="python2.7 quast.py",
                              cores=2)
        self.cores = 2
        self.to_be_removed = []
        self.copy_fasta()

    def test_make_testing_dir(self):
        """ creates temp dir for all the files created in these tests
        """
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir)
        self.assertTrue(os.path.exists(self.test_dir))

    def copy_fasta(self):
        """ make a disposable copy
        """
        copy_file(os.path.join(self.ref_dir, 'cluster1.fasta'), os.path.dirname(self.ref_fasta),
                  name="cluster1.fasta")
        self.to_be_removed.append(self.ref_fasta)

    def test_estimate_distances_smalt(self):
        """ test estimate insert disances
        """
        if os.path.exists(self.test_estimation_file):
            print("warning! existing distance esimation file!")
        est_file = estimate_distances_smalt(outfile=self.test_estimation_file,
                                            smalt_exe=self.smalt_exe,
                                            ref_genome=self.ref_fasta,
                                            fastq1=self.ref_Ffastq,
                                            fastq2=self.ref_Rfastq,
                                            cores=1,
                                            logger=logger)
        ref_smi_md5 = "a444ccbcb486a8af29736028640e87cf"  # determined manually
        ref_sma_md5 = "4ce0c8b453f2bdabd73eaf8b5ee4f376"  # determined manually
        ref_mapping_len = 9271  # mapping doesnt have exact order, so cant md5
        self.assertEqual(ref_smi_md5, md5(str(est_file + ".smi")))
        self.assertEqual(ref_sma_md5, md5(str(est_file + ".sma")))
        self.assertEqual(ref_mapping_len, file_len(est_file))
        self.to_be_removed.append(self.test_estimation_file)

    def test_ngsLib(self):
        # make a non-master object
        testlib_pe_s = ngsLib(
            name="test",
            master=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            readS0="dummy",
            ref_fasta=self.ref_fasta,
            smalt_exe=self.smalt_exe)
        self.assertFalse(os.path.exists(os.path.join(
            self.test_dir, "smalt_distance_est.sam")))
        self.assertEqual(testlib_pe_s.libtype, "pe_s")
        self.assertEqual(testlib_pe_s.readlen, None)

        # check master (ie, generate a distance file with smalt
        testlib_pe = ngsLib(
            name="test",
            master=True,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            smalt_exe=self.smalt_exe)
        self.assertTrue(os.path.exists(os.path.join(
            self.test_dir, "smalt_distance_est.sam")))
        self.assertEqual(testlib_pe.libtype, "pe")
        self.assertEqual(testlib_pe.readlen, 145.0)
        self.to_be_removed.append(testlib_pe.smalt_dist_path)

    def test_SeedGenome(self):
        gen = SeedGenome(
            max_iterations=2,
            clustered_loci_txt=self.test_loci_file,
            genbank_path=self.ref_gb,
            loci_clusters=None,
            output_root=self.test_dir)
        self.assertTrue(os.path.exists(os.path.join(self.test_dir,
                                                    "NC_011751.1.fasta")))
        self.assertTrue(os.path.exists(self.test_dir))
        self.assertTrue(isinstance(gen.seq_records[0], SeqRecord))
        with self.assertRaises(ValueError):
            SeedGenome(
                max_iterations=2,
                genbank_path=self.ref_gb,
                clustered_loci_txt=self.test_loci_file,
                loci_clusters=None,
                output_root=None)

    def test_LociMapping(self):
        with self.assertRaises(ValueError):
            LociMapping(
                name=None,
                iteration=1,
                mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        testmapping = LociMapping(
            name="test",
            iteration=1,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        self.assertTrue(os.path.isdir(testmapping.mapping_subdir))

    def test_add_coords_to_SeedGenome(self):
        gen = SeedGenome(
            max_iterations=1,
            genbank_path=self.ref_gb,
            clustered_loci_txt=self.test_loci_file,
            output_root=self.test_dir,
            logger=logger)
        gen.loci_clusters = parse_clustered_loci_file(
            filepath=gen.clustered_loci_txt,
            gb_filepath=gen.genbank_path,
            output_root=self.test_dir,
            padding=100,
            circular=False,
            logger=logger)
        add_coords_to_clusters(seedGenome=gen, logger=logger)
        print(gen.loci_clusters[0].loci_list[0].__dict__)
        self.assertEqual(gen.loci_clusters[0].loci_list[0].start_coord, 4656045)
        self.assertEqual(gen.loci_clusters[0].loci_list[0].end_coord, 4657586)

    # def test_map_to_genome_ref_smalt(self):
    #     gen = SeedGenome(
    #         max_iterations=1,
    #         genbank_path=self.ref_gb,
    #         clustered_loci_txt=self.test_loci_file,
    #         output_root=self.test_dir,
    #         logger=logger)
    #     gen.seq_ob = ngsLib(
    #         name="test",
    #         master=True,
    #         readF=self.ref_Ffastq,
    #         readR=self.ref_Rfastq,
    #         ref_fasta=self.ref_fasta,
    #         smalt_exe=self.smalt_exe)
    #     map_to_genome_ref_smalt(
    #         ngsLib=gen.seq_ob,
    #         map_results_prefix=gen.initial_map_prefix,
    #         cores=2,
    #         samtools_exe=self.samtools_exe,
    #         smalt_exe=self.smalt_exe,
    #         score_minimum=None,
    #         step=3, k=5,
    #         scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
    #         logger=logger)

    # def test_convert_bams_to_fastq_cmds(self):
    #     testmapping = LociMapping(
    #         name="test",
    #         iteration=1,
    #         mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
    #     print(testmapping.__dict__)
    #     cmds = convert_bams_to_fastq_cmds(mapping_ob=testmapping,
    #                                       ref_fasta="test_reference.fasta",
    #                                       samtools_exe=self.samtools_exe,
    #                                       which='mapped', source_ext="_bam",
    #                                       logger=logger)

    def test_generate_spades_cmds(self):
        testmapping = LociMapping(
            name="test",
            iteration=1,
            assembly_subdir=self.test_dir,
            ref_fasta=self.ref_fasta,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        testngs = ngsLib(
            name="test",
            master=True,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            smalt_exe=self.smalt_exe)

        cmd1 = generate_spades_cmd(mapping_ob=testmapping, ngs_ob=testngs,
                                   ref_as_contig='trusted',
                                   as_paired=True, addLibs="",
                                   prelim=False,
                                   k="21,33,55,77,99",
                                   spades_exe="spades.py",
                                   logger=logger)
        cmd1_ref = str("spades.py --careful -k 21,33,55,77,99 --pe1-1 {0} " +
                       "--pe1-2 {1}  --trusted-contigs {2}  -o {3}").format(
                           self.ref_Ffastq, self.ref_Rfastq, self.ref_fasta,
                           self.output_root)
        self.assertEqual(cmd1, cmd1_ref)
        self.to_be_removed.append(testngs.smalt_dist_path)

    def test_run_final_assembliies(self):
        gen = SeedGenome(
            max_iterations=1,
            genbank_path=self.ref_gb,
            clustered_loci_txt=self.test_loci_file,
            output_root=self.test_dir,
            logger=logger)

        final_cmds = run_final_assemblies(args=self.args,
                                          seedGenome=gen, logger=logger)
        for i in final_cmds:
            print(i)
        pass


    # def test_convert_spades_cmds(self):
    #     get_extract_convert_spades_cmds(mapping_ob, fetch_mates, samtools_exe,
    #                                     spades_exe, ref_as_contig, logger)
    # def test_partition_mapped_reads(self):
    #     gen = SeedGenome(
    #         max_iterations=1,
    #         genbank_path=self.ref_gb,
    #         clustered_loci_txt=self.test_loci_file,
    #         output_root=self.test_dir,
    #         # initial_map_prefix=os.path.join(self.test_dir, "LociMapping"),
    #         logger=logger)
    #     gen.ngs_ob = ngsLib(
    #         name="test",
    #         master=True,
    #         readF=self.ref_Ffastq,
    #         readR=self.ref_Rfastq,
    #         readS0=None,
    #         ref_fasta=gen.ref_fasta,
    #         smalt_dist_path=None,
    #         # readlen=None,
    #         smalt_exe=self.smalt_exe)
    #     gen.loci_clusters = parse_clustered_loci_file(
    #         filepath=gen.riboSelect_path,
    #         output_root=self.test_dir,
    #         gb_filepath=gen.genbank_path,
    #         padding=100,
    #         circular=False,
    #         logger=logger)
    #     add_coords_to_clusters(seedGenome=gen, logger=logger)
    #     map_to_genome_ref_smalt(
    #         ref=gen.ref_fasta,
    #         ngsLib=gen.ngs_ob,
    #         map_results_prefix=gen.initial_map_prefix,
    #         cores=2,
    #         samtools_exe=self.samtools_exe,
    #         smalt_exe=self.smalt_exe,
    #         score_minimum=None,
    #         step=3, k=5,
    #         scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
    #         logger=logger)
    #     # add output root to each
    #     partition_mapping(
    #         seedGenome=gen,
    #         samtools_exe=self.samtools_exe,
    #         flank=[0, 0],
    #         logger=logger)
    #     logger.warning("running without multiprocessing!")
    #     for cluster in gen.loci_clusters:
    #         assemble_iterative_mapping(
    #             clu=cluster,
    #             master_ngs_ob=gen.ngs_ob,
    #             args=self,
    #             nseqs=len(gen.loci_clusters),
    #             fetch_mates=False,
    #             include_short_contigs=True,
    #             max_iterations=2,
    #             min_contig_len=3000,
    #             target_len=7000,
    #             samtools_exe=self.samtools_exe,
    #             keep_unmapped_reads=False,
    #             logger=logger)

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
