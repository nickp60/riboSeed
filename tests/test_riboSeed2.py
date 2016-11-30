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
# import subprocess
import os
import unittest
# import multiprocessing

# from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import Namespace

# I hate this line but it works :(
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))


from pyutilsnrw.utils3_5 import md5, file_len, copy_file, get_number_mapped

from riboSeed.riboSeed2 import SeedGenome, ngsLib,  LociMapping, \
    map_to_genome_ref_smalt, add_coords_to_clusters, partition_mapping, \
    convert_bams_to_fastq_cmds, check_smalt_full_install,\
    generate_spades_cmd, estimate_distances_smalt, run_final_assemblies,\
    check_libs_before_mapping, make_faux_genome

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
        self.spades_exe = "spades.py"
        self.quast_exe = "quast.py"
        self.quast_python_exe = "python2"
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
        self.maxDiff = 2000
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
        testlib_s = ngsLib(
            name="test",
            master=False,
            readF=None,
            readR=None,
            readS0=self.ref_Ffastq,
            ref_fasta=self.ref_fasta,
            smalt_exe=self.smalt_exe)
        self.assertEqual(testlib_s.libtype, "s_1")
        self.assertEqual(testlib_s.readlen, 145.0)
        # test unnamed fails
        with self.assertRaises(ValueError):
            ngsLib(
                name=None,
                master=False,
                readF=self.ref_Ffastq,
                readR=self.ref_Rfastq,
                readS0="dummy",
                ref_fasta=self.ref_fasta,
                smalt_exe=self.smalt_exe)
        # self.assertFalse(os.path.exists(os.path.join(
        #     self.test_dir, "smalt_distance_est.sam")))
        self.assertEqual(testlib_pe_s.libtype, "pe_s")
        self.assertEqual(testlib_pe_s.readlen, None)
        # test fails with singe PE file
        with self.assertRaises(ValueError):
            ngsLib(
                name=None,
                master=False,
                readF=self.ref_Ffastq,
                readR=None,
                readS0="dummy",
                ref_fasta=self.ref_fasta,
                smalt_exe=self.smalt_exe)

        # check master (ie, generate a distance file with smalt
        testlib_pe = ngsLib(
            name="test",
            master=True,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            smalt_exe=self.smalt_exe,
            logger=logger)
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

    def test_check_smalt_full_install(self):
        """ TODO: how would I test this?
        """
        smalttestdir = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                                    "sample_data",
                                    "smalt_test", "")
        ref = os.path.join(smalttestdir, "ref_to_test_bambamc.fasta")
        index = os.path.join(smalttestdir, "test_index")
        test_bam = os.path.join(smalttestdir, "test_mapping.bam")
        test_reads = os.path.join(smalttestdir, "reads_to_test_bambamc.fastq")

        check_cmds = check_smalt_full_install(self.smalt_exe, logger=logger)
        check_cmds_ref = [
            "{0} index {1} {2}".format(
                self.smalt_exe, index, ref),
            "{0} map -f bam -o {1} {2} {3}".format(
                self.smalt_exe, test_bam, index, test_reads)]
        for index, cmd in enumerate(check_cmds):
            self.assertEqual(cmd, check_cmds_ref[index])

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
        self.assertEqual(
            gen.loci_clusters[0].loci_list[0].start_coord, 4656045)
        self.assertEqual(
            gen.loci_clusters[0].loci_list[0].end_coord, 4657586)

    def test_convert_bams_to_fastq_cmds(self):
        testmapping = LociMapping(
            name="test",
            iteration=1,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        cmd, ngs = convert_bams_to_fastq_cmds(mapping_ob=testmapping,
                                              ref_fasta="test_reference.fasta",
                                              samtools_exe=self.samtools_exe,
                                              which='mapped',
                                              source_ext="_bam",
                                              logger=logger)
        cmd_ref = "{0} fastq {1} -1 {2} -2 {3} -s {4}".format(
            self.samtools_exe, testmapping.mapped_bam,
            os.path.join(testmapping.mapping_subdir,
                         "test_iteration_1_mappedreadF.fastq"),
            os.path.join(testmapping.mapping_subdir,
                         "test_iteration_1_mappedreadR.fastq"),
            os.path.join(testmapping.mapping_subdir,
                         "test_iteration_1_mappedreadS.fastq"), )
        self.assertEqual(cmd[0], cmd_ref)

    def test_generate_spades_cmds(self):
        """Question: why the heck is there a space before the -o?
        I troubleshot this for two hours to no avail.  Next!
        """
        testmapping = LociMapping(
            name="test",
            iteration=1,
            assembly_subdir=self.test_dir,
            ref_fasta=self.ref_fasta,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        testngs = ngsLib(
            name="test",
            master=False,
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
        cmd1_ref = str(
            "spades.py --careful -k 21,33,55,77,99 --pe1-1 {0} " +
            "--pe1-2 {1} --trusted-contigs {2}  -o {3}").format(
            self.ref_Ffastq, self.ref_Rfastq, self.ref_fasta, self.test_dir)
        self.assertEqual(cmd1, cmd1_ref)

    def test_lib_check(self):
        empty_file = os.path.join(self.test_dir, "test_not_real_file")
        # make an empty file
        with open(empty_file, 'w') as ef:
            pass
        ngs_ob = ngsLib(
            name="test",
            master=False,
            readF=self.ref_Ffastq,
            readR=empty_file,
            ref_fasta=self.ref_fasta,
            smalt_exe=self.smalt_exe)
        check_libs_before_mapping(ngsLib=ngs_ob, logger=logger)
        self.assertTrue(ngs_ob.readR is None)
        self.to_be_removed.append(empty_file)

    def test_map_with_smalt(self):
        # becasue multiple mapping are assingmed randomly (pseudorandomly?),
        # this accepts a range of expected results
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

        map_to_genome_ref_smalt(mapping_ob=testmapping, ngsLib=testngs,
                                cores=4, samtools_exe=self.samtools_exe,
                                smalt_exe=self.smalt_exe, score_minimum=45,
                                scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
                                step=3, k=5, logger=logger)
        mapped_str = get_number_mapped(testmapping.pe_map_bam,
                                       samtools_exe=self.samtools_exe)
        nmapped = int(mapped_str[0:5])
        nperc = float(mapped_str[-13:-8])
        print(nmapped)
        print(nperc)
        self.assertTrue(12575 < nmapped < 12635)
        self.assertTrue(34.6 < nperc < 34.8)

    def test_make_faux_genome(self):
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

        for loci in gen.loci_clusters:
            loci.keep_contig = True
            loci.mappings[0] = LociMapping(
                name="test",
                iteration=1,
                assembly_subdir=self.test_dir,
                ref_fasta=self.ref_fasta,
                assembled_contig=self.ref_fasta,
                mapping_subdir=os.path.join(self.test_dir, "LociMapping_for_text_faux_genome"))

        # map_to_genome_ref_smalt(
        #     mapping_ob=gen.iter_mapping_list[0],
        #     ngsLib=gen.master_ngs_ob,
        #     cores=4,
        #     samtools_exe=self.samtools_exe,
        #     smalt_exe=self.smalt_exe,
        #     score_minimum=None,
        #     step=3, k=5,
        #     scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
        #     logger=logger)
        # partition_mapping(seedGenome=gen,
        #                   logger=logger,
        #                   samtools_exe=self.samtools_exe,
        #                   flank="50:50",
        #                   cluster_list=gen.loci_clusters)

        make_faux_genome(cluster_list=gen.loci_clusters, seedGenome=gen,
                         iteration=1, output_root=self.test_dir, nbuff=10000,
                         logger=logger)

    def test_run_final_assembliies(self):
        gen = SeedGenome(
            max_iterations=1,
            genbank_path=self.ref_gb,
            clustered_loci_txt=self.test_loci_file,
            output_root=self.test_dir,
            assembled_seeds="tralalalala",
            logger=logger,
            master_ngs_ob=ngsLib(
                name="test",
                master=False,
                readF=self.ref_Ffastq,
                readR=self.ref_Rfastq,
                ref_fasta=self.ref_fasta,
                smalt_exe=self.smalt_exe))
        gen.ref_fasta = self.ref_fasta
        final_cmds, quast_reports = run_final_assemblies(
            seedGenome=gen, spades_exe=self.spades_exe,
            quast_exe=self.quast_exe, quast_python_exe=self.quast_python_exe,
            skip_control=False, kmers="33,77,99", logger=logger)
        final_cmds_ref = [
            str(
                "{0} --careful -k 33,77,99 --pe1-1 {1} " +
                "--pe1-2 {2} --trusted-contigs {3}  -o {4}"
            ).format(
                self.spades_exe, self.ref_Ffastq, self.ref_Rfastq,
                gen.assembled_seeds,
                os.path.join(self.test_dir, "final_de_fere_novo_assembly")),
            str(
                '{0} {1} tralalalala -R {2} -o {3}'
            ).format(
                self.quast_python_exe, self.quast_exe,
                self.ref_fasta,
                os.path.join(self.test_dir, "quast_de_fere_novo")),
            str(
                "{0} --careful -k 33,77,99 --pe1-1 {1} " +
                "--pe1-2 {2}   -o {3}"
            ).format(
                self.spades_exe, self.ref_Ffastq, self.ref_Rfastq,
                os.path.join(self.test_dir, "final_de_novo_assembly")),
            str(
                '{0} {1} tralalalala -R {2} -o {3}'
            ).format(
                self.quast_python_exe, self.quast_exe,
                self.ref_fasta,
                os.path.join(self.test_dir, "quast_de_novo"))]
        for i, ref in enumerate(final_cmds_ref):
            self.assertEqual(final_cmds[i], ref)

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
