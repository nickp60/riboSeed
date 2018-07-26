# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import shutil
from unittest.mock import patch, create_autospec
import os
import unittest
import time
import subprocess

from Bio import SeqIO
from argparse import Namespace


from riboSeed.shared_methods import md5, get_number_mapped, check_version_from_cmd

from riboSeed.classes import SeedGenome, NgsLib, LociMapping, Exes
from riboSeed.riboSeed import \
    map_to_genome_ref_smalt, map_to_genome_ref_bwa, \
    add_coords_to_clusters, partition_mapping, \
    convert_bam_to_fastqs_cmd, get_smalt_full_install_cmds,\
    generate_spades_cmd, get_final_assemblies_cmds,\
    nonify_empty_lib_files, make_faux_genome, \
    evaluate_spades_success, prepare_next_mapping, make_mapped_partition_cmds,\
    make_unmapped_partition_cmds, make_quick_quast_table, \
    check_fastqs_len_equal, parse_subassembly_return_code, \
    make_spades_empty_check, get_samtools_depths, \
    decide_proceed_to_target, get_rec_from_generator, \
    check_kmer_vs_reads, make_samtools_depth_cmds, \
    parse_samtools_depth_results, make_modest_spades_cmd, get_bam_AS, \
    pysam_extract_reads, define_score_minimum, bool_run_quast, \
    make_quast_command, check_genbank_for_fasta, \
    filter_bam_AS, make_bwa_map_cmds, set_ref_as_contig, \
    make_get_consensus_cmds, exclude_subassembly_based_on_coverage, \
    convert_fastq_to_fasta, get_fasta_consensus_from_BAM, \
    fiddle_with_spades_exe, check_spades_extra_library_input

from riboSeed.shared_methods import parse_clustered_loci_file

sys.dont_write_bytecode = True
logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class riboSeedShallow(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_riboSeed_tests")
        self.spades_dir = os.path.join(os.path.dirname(__file__),
                                       "output_riboSeed_tests",
                                       "SPAdes_results")
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.ref_gb = os.path.join(self.ref_dir,
                                   'NC_011751.1.gb')
        self.ref_tiny_gb = os.path.join(self.ref_dir,
                                        'scannedScaffolds.gb')
        self.ref_fasta = os.path.join(self.test_dir,
                                      'cluster1.fasta')
        self.good_contig = os.path.join(self.ref_dir,
                                        'contigs.fasta')
        self.short_contig = os.path.join(self.ref_dir,
                                         'contigs.fasta')
        self.ref_Ffastq = os.path.join(self.ref_dir,
                                       'toy_reads1.fq')
        self.ref_Rfastq = os.path.join(self.ref_dir,
                                       'toy_reads2.fq')
        self.ref_Sfastq = os.path.join(self.ref_dir,
                                       'toy_readsS.fq')
        self.ref_bam_prefix = os.path.join(self.ref_dir,
                                           'test_bam_to_fastq')
        self.ref_fastq_with_iupac_bases = os.path.join(
            self.ref_dir,
            "riboSeed_references",
            "file_with_iupac_bases.fastq")
        self.smalt_exe = "smalt"
        self.bwa_exe = "bwa"
        self.samtools_exe = "samtools"
        self.spades_exe = "spades.py"
        self.quast_exe = "quast"
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
                              spades_exe="spades.py",
                              cores=2)
        self.cores = 2
        self.maxDiff = 2000
        self.testmapping =  LociMapping(
            name="test",
            iteration=1,
            assembly_subdir=self.test_dir,
            ref_fasta=self.ref_fasta,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        self.testngs1 = NgsLib(
            name="test",
            master=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)
        self.testngs2 = NgsLib(
            name="test",
            master=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            readS0=self.ref_Sfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)

        self.to_be_removed = []
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)
        self.copy_fasta()
        self.startTime = time.time()
        # self.print_headsup()

    def copy_fasta(self):
        """ make a disposable copy of the cluster fasta
        """
        shutil.copy(os.path.join(self.ref_dir, 'cluster1.fasta'),
                    self.ref_fasta)
        self.to_be_removed.append(self.ref_fasta)

    def print_headsup(self):
        """ Warn about some expected error messages
        """
        print("""
        Thanks for running these tests!! You will see several logger warning
        messages that you should IGNORE. **sigh of relief** They verify error
        reporting in the actual function, and cannot be removed. You should
        see the following errors if all went properly:
        \t ERROR:root:--target_len is set to invalid value! Must be a decimal greater than zero, ie where 1.1 would be 110% of the original sequence length.
        \t ERROR:root:--target_len is set to invalid value! Must be a decimal greater than zero, ie where 1.1 would be 110% of the original sequence length.
        \t ERROR:root:We dont reccommend seeding to lengths greater than5x original seed length. Try between 0.5 and 1.5.  If you are setting a target number of bases, it  must be greater than 50
        \t WARNING:root:The first iteration's assembly's best contig is not greater than length set by --min_assembly_len. Assembly will likely fail if the contig does not meet the length of the seed
        \t WARNING:root:Continuing, but if this occurs for more than one seed, we reccommend  you abort and retry with longer seeds, a different ref, or re-examine the riboSnag clustering
        \t WARNING:root:The first iteration's assembly's best contig is not greater than length set by --min_assembly_len. Assembly will likely fail if the contig does not meet the length of the seed
        \t WARNING:root:SEED_cluster-7-iter1: No output from SPAdes this time around...
        \t WARNING:root:read file readR is empty and will not be used for mapping!
        """)

    def test_make_spades_empty_check(self):
        """ construct spades command that check for file presense
        this is useful when multiprocessing and unable to check before
        sending the command out
        """
        lib1 = "/path/to/lib1.fastq"
        lib2 = "/path/to/lib2.fastq"
        lib3 = "/path/to/lib3.fastq"
        cmd = "spades command"
        fullcmd = make_spades_empty_check(
            liblist=[lib1, lib2, lib3],
            cmd=cmd,
            logger=logger)
        self.assertEqual(
            str("if [ -s {0} ] && [ -s {1} ] && [ -s {2} ] ; " +
                "then {3} ; else echo 'input lib not found, " +
                "skipping this SPAdes call' ; " +
                "fi").format(lib1, lib2, lib3, cmd),
            fullcmd)

    def test_LociMapping(self):
        """ Does LociMapping instatiate correctly
        """
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

    def test_check_kmer_bad_input(self):
        """ does kmer checker fail with a bad input
        """
        with self.assertRaises(ValueError):
            # non integer
            check_kmer_vs_reads(k="33,22,55,spinach,77",
                                readlen=50, min_diff=2, logger=logger)
        with self.assertRaises(ValueError):
            # non comma delim
            check_kmer_vs_reads(k="33,22;55,77",
                                readlen=50, min_diff=2, logger=logger)

    def test_check_kmer_skip_auto(self):
        """ does kmer checker ignore 'auto' properly
        """
        self.assertEqual("auto",
                         check_kmer_vs_reads(k="auto",
                                             readlen=50,
                                             min_diff=2, logger=logger))

    def test_check_kmer_filter_big(self):
        """ does kmer checker remove too big ks
        """
        self.assertEqual("21,33",
                         check_kmer_vs_reads(k="21,33,55,77",
                                             readlen=50, min_diff=2,
                                             logger=logger))

    def test_check_kmer_filter_close(self):
        """ does kmer checker remove too close k
        """
        self.assertEqual("21,33,55",
                         check_kmer_vs_reads(k="21,33,55,77",
                                             readlen=79, min_diff=2,
                                             logger=logger))

    def test_check_kmer_filter_even(self):
        """ does kmer checker remove even k
        """
        self.assertEqual("21,33,55",
                         check_kmer_vs_reads(k="21,33,55,77,90",
                                             readlen=79, min_diff=2,
                                             logger=logger))

    def test_check_kmer_good(self):
        """ does kmer checker work when all is well
        """
        self.assertEqual("21,33,55,77,127",
                         check_kmer_vs_reads(k="21,33,55,77,127",
                                             readlen=250, min_diff=2,
                                             logger=logger))

    def test_get_smalt_full_install_cmds(self):
        """ can we make the commands to check smalt instaltion
        """
        smalttestdir = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                                    "sample_data",
                                    "smalt_test", "")
        ref = os.path.join(smalttestdir, "ref_to_test_bambamc.fasta")
        index = os.path.join(smalttestdir, "test_index")
        test_bam = os.path.join(smalttestdir, "test_mapping.bam")
        test_reads = os.path.join(smalttestdir, "reads_to_test_bambamc.fastq")

        check_cmds = get_smalt_full_install_cmds(self.smalt_exe, logger=logger)
        check_cmds_ref = [
            "{0} index {1} {2}".format(
                self.smalt_exe, index, ref),
            "{0} map -f bam -o {1} {2} {3}".format(
                self.smalt_exe, test_bam, index, test_reads)]
        for index, cmd in enumerate(check_cmds):
            self.assertEqual(cmd, check_cmds_ref[index])

    def test_get_rec_from_generator(self):
        """ can we get a record of interest from a SeqIO generator
        tests get_rec_from_generator, which replaces
        get_genbank_rec_from_multigb
        """
        recsgen = None
        with open(self.ref_gb, "r") as gbf:
            global reclist
            reclist = list(SeqIO.parse(gbf, "genbank"))
            return reclist

        def refreshGen():
            recsgen = SeqIO.parse(self.ref_gb, "genbank")

        refreshGen()
        record = get_rec_from_generator(recordID='NC_011751.1',
                                        method=refreshGen,
                                        gen=recsgen)
        # this check that it got refreshed
        self.assertEqual(record.id, tuple(recsgen)[0].id)
        self.assertEqual(record.id, reclist[0].id)
        with self.assertRaises(ValueError):
            get_rec_from_generator(recordID='NC_011751.X',
                                   method=refreshGen,
                                   gen=recsgen)

    def test_convert_bam_to_fastqs_cmd(self):
        """ can we construct proper samtools fast cmds
        """
        testmapping = LociMapping(
            name="test",
            iteration=1,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        cmd, ngs = convert_bam_to_fastqs_cmd(mapping_ob=testmapping,
                                             ref_fasta="test_reference.fasta",
                                             samtools_exe=self.samtools_exe,
                                             which='mapped',
                                             source_ext="_bam",
                                             logger=logger)
        cmd_ref = "{0} fastq {1} -1 {2} -2 {3} -s {4}".format(
            self.samtools_exe, testmapping.mapped_bam,
            os.path.join(testmapping.mapping_subdir,
                         "test_mappedreadF.fastq"),
            os.path.join(testmapping.mapping_subdir,
                         "test_mappedreadR.fastq"),
            os.path.join(testmapping.mapping_subdir,
                         "test_mappedreadS.fastq"), )
        self.assertEqual(cmd, cmd_ref)

    def test_generate_spades_cmds_pe_prelim(self):
        """ PE reads, prelim
        """
        cmd1 = generate_spades_cmd(mapping_ob=self.testmapping,
                                   ngs_ob=self.testngs1,
                                   ref_as_contig='trusted',
                                   as_paired=True, addLibs="",
                                   prelim=True,
                                   k="21,33,55,77,99",
                                   python_exe="python3.5",
                                   spades_exe="spades.py",
                                   logger=logger)
        cmd1_ref = str(
            "python3.5 spades.py --only-assembler --cov-cutoff off --sc --careful " +
            "-k 21,33,55,77,99 --pe1-1 {0} " +
            "--pe1-2 {1} --trusted-contigs {2}  -o {3}"
        ).format(
            self.ref_Ffastq, self.ref_Rfastq, self.ref_fasta, self.test_dir)
        self.assertEqual(cmd1, cmd1_ref)

    def test_generate_spades_cmds_pe_s(self):
        """ PE with singletons
        """
        cmd2_ref = str(
            "python3.5 spades.py --careful -k 21,33,55,77,99 --pe1-1 {0} " +
            "--pe1-2 {1} --pe1-s {2} --trusted-contigs {3}  -o {4}"
        ).format(
            self.ref_Ffastq, self.ref_Rfastq, self.ref_Sfastq,
            self.ref_fasta, self.test_dir)
        cmd2 = generate_spades_cmd(mapping_ob=self.testmapping,
                                   ngs_ob=self.testngs2,
                                   ref_as_contig='trusted',
                                   as_paired=True, addLibs="",
                                   prelim=False,
                                   k="21,33,55,77,99",
                                   python_exe="python3.5",
                                   spades_exe="spades.py",
                                   logger=logger)
        self.assertEqual(cmd2, cmd2_ref)

    def test_generate_spades_cmds_pe_as_s(self):
        """ PE as sngle libraries
        """
        cmd3 = generate_spades_cmd(mapping_ob=self.testmapping,
                                   ngs_ob=self.testngs1,
                                   ref_as_contig='trusted',
                                   as_paired=False, addLibs="",
                                   prelim=False,
                                   k="21,33,55,77,99",
                                   python_exe="python3.5",
                                   spades_exe="spades.py",
                                   logger=logger)
        cmd3_ref = str(
            "python3.5 spades.py --careful -k 21,33,55,77,99 --pe1-s {0} " +
            "--pe2-s {1} --trusted-contigs {2}  -o {3}"
        ).format(
            self.ref_Ffastq, self.ref_Rfastq, self.ref_fasta, self.test_dir)
        self.assertEqual(cmd3, cmd3_ref)

    def test_generate_spades_cmds_pe_and_s_as_s(self):
        """ PE and singletons all as single libraries
        """
        cmd4 = generate_spades_cmd(mapping_ob=self.testmapping,
                                   ngs_ob=self.testngs2,
                                   ref_as_contig='trusted',
                                   as_paired=False, addLibs="",
                                   prelim=False,
                                   k="21,33,55,77,99",
                                   python_exe="python3.5",
                                   spades_exe="spades.py",
                                   logger=logger)
        cmd4_ref = str(
            "python3.5 spades.py --careful -k 21,33,55,77,99 --pe1-s {0} " +
            "--pe2-s {1} --pe3-s {4}  --trusted-contigs {2}  -o {3}"
        ).format(
            self.ref_Ffastq, self.ref_Rfastq, self.ref_fasta, self.test_dir,
        self.ref_Sfastq)
        self.assertEqual(cmd4, cmd4_ref)

    def test_parse_subassembly_return_code_0(self):
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
        gen.loci_clusters[0].assembly_success = 0
        parse_subassembly_return_code(cluster=gen.loci_clusters[0],
                                      logger=logger)
        self.assertTrue(gen.loci_clusters[0].continue_iterating)
        self.assertTrue(gen.loci_clusters[0].keep_contigs)

    def test_parse_subassembly_return_code_1(self):
        """ this is depreciated, and should raise a valueError
        """
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
        gen.loci_clusters[0].assembly_success = 1
        with self.assertRaises(ValueError):
            parse_subassembly_return_code(cluster=gen.loci_clusters[0],
                                          logger=logger)

    def test_parse_subassembly_return_code_2(self):
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
        gen.loci_clusters[0].assembly_success = 2
        parse_subassembly_return_code(cluster=gen.loci_clusters[0],
                                      logger=logger)
        self.assertFalse(gen.loci_clusters[0].continue_iterating)
        self.assertFalse(gen.loci_clusters[0].keep_contigs)

    def test_parse_subassembly_return_code_3(self):
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
        gen.loci_clusters[0].assembly_success = 3
        parse_subassembly_return_code(cluster=gen.loci_clusters[0],
                                      logger=logger)
        self.assertFalse(gen.loci_clusters[0].continue_iterating)
        self.assertFalse(gen.loci_clusters[0].keep_contigs)

    def test_check_fastqs_len_equal(self):
        failed = False
        try:
            check_fastqs_len_equal(self.ref_fasta, self.ref_fasta)
        except:
            failed = True
        self.assertFalse(failed)

    def test_check_fastqs_len_equal_fail(self):
        with self.assertRaises(ValueError):
            check_fastqs_len_equal(self.ref_fasta, self.ref_gb)

    def test_evaluate_spades(self):
        """ Test all possible combinations of spades return codes
        """
        shutil.copyfile(self.good_contig,
                        os.path.join(self.test_dir,
                                     "contigs.fasta"))
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
        testmapping = LociMapping(
            name="test",
            iteration=1,
            assembly_subdir=self.test_dir,
            ref_fasta=self.ref_fasta,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        testmapping2 = LociMapping(
            name="test",
            iteration=1,
            assembly_subdir_needed=False,
            assembly_subdir="nonexistantdir",
            ref_fasta=self.ref_fasta,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        clu = gen.loci_clusters[0]
        clu.mappings.append(testmapping)
        code_0 = evaluate_spades_success(
            clu,
            mapping_ob=testmapping,
            proceed_to_target=False, target_len=None,
            min_assembly_len=5000,
            include_short_contigs=False, keep_best_contig=True,
            seqname='', logger=logger)
        self.assertEqual(code_0, 0)
        # code_1 = evaluate_spades_success(
        #     clu,
        #     mapping_ob=testmapping,
        #     proceed_to_target=False, target_len=None,
        #     read_len=300,
        #     min_assembly_len=10000,
        #     include_short_contigs=True, keep_best_contig=True,
        #     seqname='', logger=logger)
        # self.assertEqual(code_1, 1)
        # code_1a = evaluate_spades_success(  # test proceed_to_target
        #     clu,
        #     mapping_ob=testmapping,
        #     proceed_to_target=True, target_len=6000,
        #     read_len=300,
        #     min_assembly_len=5000,
        #     include_short_contigs=True, keep_best_contig=True,
        #     seqname='', logger=logger)
        # self.assertEqual(code_1a, 1)
        code_2 = evaluate_spades_success(
            clu,
            mapping_ob=testmapping,
            proceed_to_target=False, target_len=None,
            min_assembly_len=10000,
            include_short_contigs=False, keep_best_contig=True,
            seqname='', logger=logger)
        self.assertEqual(code_2, 2)
        code_3 = evaluate_spades_success(
            clu,
            mapping_ob=testmapping2,
            proceed_to_target=False, target_len=None,
            min_assembly_len=5000,
            include_short_contigs=False, keep_best_contig=True,
            seqname='', logger=logger)
        self.assertEqual(code_3, 3)

    def test_make_quast_table(self):
        paths = [os.path.join(self.ref_dir, "quast_1", "report.tsv"),
                 os.path.join(self.ref_dir, "quast_2", "report.tsv")]
        quast_dict = make_quick_quast_table(
            pathlist=paths, write=False, writedir=None, logger=logger)
        first10_printout = [
            "# N's per 100 kbp: ['0.00', '0.00']",
            "# contigs: ['2', '11']",
            "# contigs (>= 0 bp): ['3', '20']",
            "# contigs (>= 1000 bp): ['2', '10']",
            "# contigs (>= 10000 bp): ['2', '3']",
            "# contigs (>= 25000 bp): ['1', '0']",
            "# contigs (>= 5000 bp): ['2', '7']",
            "# contigs (>= 50000 bp): ['1', '0']",
            "# indels per 100 kbp: ['3.98', '0.00']",
            "# local misassemblies: ['0', '0']",
            "# misassembled contigs: ['0', '0']"]
        counter = 0
        for k, v in sorted(quast_dict.items()):
            self.assertEqual(first10_printout[counter],
                             "{0}: {1}".format(k, v))
            if counter == 10:
                break
            else:
                counter = counter + 1

    def test_partition_mapping(self):
        # this mostly does system calls; cant really test smoothly
        pass

    def test_mapped_partition_cmds(self):
        gen = SeedGenome(
            max_iterations=1,
            genbank_path=self.ref_gb,
            clustered_loci_txt=self.test_loci_file,
            output_root=self.test_dir,
            logger=logger)
        testmapping = LociMapping(
            name="test",
            iteration=1,
            assembly_subdir=self.test_dir,
            ref_fasta=self.ref_fasta,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        gen.loci_clusters = parse_clustered_loci_file(
            filepath=gen.clustered_loci_txt,
            gb_filepath=gen.genbank_path,
            output_root=self.test_dir,
            padding=1000,
            circular=False,
            logger=logger)
        # add_coords_to_clusters(seedGenome=gen, logger=logger)
        gen.loci_clusters[0].global_start_coord = 5000
        gen.loci_clusters[0].global_end_coord = 10000
        # print(gen.loci_clusters[0].__dict__)
        cmds, region = make_mapped_partition_cmds(
            cluster=gen.loci_clusters[0],
            mapping_ob=testmapping,
            seedGenome=gen,
            samtools_exe=self.samtools_exe)
        ref_cmds = [
            '{0} sort {1} > {2}'.format(
                self.samtools_exe,
                os.path.join(
                    self.test_dir,
                    "NC_011751.1_mapping_for_iteration_0",
                    "NC_011751.1_mapping_iteration_0.bam"),
                os.path.join(
                    self.test_dir,
                    "NC_011751.1_mapping_for_iteration_0",
                    "NC_011751.1_mapping_iteration_0_sorted.bam")),
            "{0} index {1}".format(
                self.samtools_exe,
                os.path.join(
                    self.test_dir,
                    "NC_011751.1_mapping_for_iteration_0",
                    "NC_011751.1_mapping_iteration_0_sorted.bam")),
            '{0} view -o {1} {2} NC_011751.1:5000-10000'.format(
                self.samtools_exe,
                os.path.join(self.test_dir,
                             "LociMapping", "test.bam"),
                os.path.join(
                    self.test_dir,
                    "NC_011751.1_mapping_for_iteration_0",
                    "NC_011751.1_mapping_iteration_0_sorted.bam"))]
        # print(region)
        self.assertEqual(region, "NC_011751.1:5000-10000")
        for i, cmd in enumerate(cmds):
            self.assertEqual(ref_cmds[i], cmd)

    def test_unmapped_partition_cmds(self):
        gen = SeedGenome(
            max_iterations=1,
            genbank_path=self.ref_gb,
            clustered_loci_txt=self.test_loci_file,
            output_root=self.test_dir,
            logger=logger)

        unmapped_cmds = make_unmapped_partition_cmds(
            mapped_regions=["test:1-6", "test:3-77"],
            samtools_exe=self.samtools_exe, seedGenome=gen)
        mapped_txt = os.path.join(
            self.test_dir,
            "NC_011751.1_mapping_for_iteration_0",
            "NC_011751.1_mapping_iteration_0_mapped.txt")

        ref_unmapped_cmds = [
            "{0} view -o {1} -h {2}".format(
                self.samtools_exe,
                os.path.join(self.test_dir,
                             "NC_011751.1_mapping_for_iteration_0",
                             "NC_011751.1_mapping_iteration_0.sam"),
                os.path.join(
                    self.test_dir,
                    "NC_011751.1_mapping_for_iteration_0",
                    "NC_011751.1_mapping_iteration_0.bam")),
            "{0} view {1} test:1-6 | cut -f1 >> {2}".format(
                self.samtools_exe,
                os.path.join(
                    self.test_dir,
                    "NC_011751.1_mapping_for_iteration_0",
                    "NC_011751.1_mapping_iteration_0_sorted.bam"),
                mapped_txt),
            "{0} view {1} test:3-77 | cut -f1 >> {2}".format(
                self.samtools_exe,
                os.path.join(
                    self.test_dir,
                    "NC_011751.1_mapping_for_iteration_0",
                    "NC_011751.1_mapping_iteration_0_sorted.bam"),
                mapped_txt),
            "sort -u {0} -o {0}".format(mapped_txt)]
        for i, ref in enumerate(ref_unmapped_cmds):
            self.assertEqual(unmapped_cmds[i], ref)
        # self.to_be_removed.append(mapped_txt)

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
        for clu in gen.loci_clusters:
            clu.keep_contig = True
            clu.mappings.append(LociMapping(
                name="test",
                iteration=1,
                assembly_subdir=self.test_dir,
                ref_fasta=self.ref_fasta,
                assembled_contig=self.ref_fasta,
                mapping_subdir=os.path.join(
                    self.test_dir, "LociMapping_for_text_faux_genome")))

        make_faux_genome(cluster_list=gen.loci_clusters, seedGenome=gen,
                         iteration=1, output_root=self.test_dir, nbuff=10000,
                         logger=logger)
        # add a test here

    def test_get_final_assemblies_cmds(self):
        gen = SeedGenome(
            max_iterations=1,
            genbank_path=self.ref_gb,
            clustered_loci_txt=self.test_loci_file,
            output_root=self.test_dir,
            assembled_seeds="tralalalala",
            logger=logger,
            master_ngs_ob=NgsLib(
                name="test",
                master=False,
                readF=self.ref_Ffastq,
                readR=self.ref_Rfastq,
                ref_fasta=self.ref_fasta,
                mapper_exe=self.smalt_exe))
        test_exes = Exes(python="/bin/python3.5",
                         samtools=self.samtools_exe,
                         quast=self.quast_exe,
                         smalt=self.smalt_exe,
                         spades=self.spades_exe,
                         bwa=self.bwa_exe,
                         check=False,
                         method="bwa")
        # I want this  generalized, so replace the acual path with spades.py
        test_exes.spades = "spades.py"
        test_exes.quast = "quast"
        gen.ref_fasta = self.ref_fasta
        with patch.object(sys, 'version_info') as v_info:
            v_info.minor = 5
            final_cmds, quast_reports = \
                get_final_assemblies_cmds(
                    ref_as_contig="trusted",
                    cores=4,
                    serialize=True,
                    additional_libs="-s1 path/to/single_reads.fastq",
                    memory=8,
                    seedGenome=gen, exes=test_exes,
                    skip_control=False, kmers="33,77,99", logger=logger)
        final_spades_cmds_ref = [
            str(
                "/bin/python3.5 {0} -t 4 -m 8 --careful -k 33,77,99 --pe1-1 {1} " +
                "--pe1-2 {2} --trusted-contigs {3}  -o {4} -s1 path/to/single_reads.fastq"
            ).format(
                self.spades_exe, self.ref_Ffastq, self.ref_Rfastq,
                gen.assembled_seeds,
                os.path.join(self.test_dir, "final_de_fere_novo_assembly")),
            str(
                "/bin/python3.5 {0} -t 4 -m 8 --careful -k 33,77,99 --pe1-1 {1} " +
                "--pe1-2 {2}   -o {3} -s1 path/to/single_reads.fastq"
            ).format(
                self.spades_exe, self.ref_Ffastq, self.ref_Rfastq,
                os.path.join(self.test_dir, "final_de_novo_assembly"))]
        final_quast_cmds_ref = [
            str(
                '{0} {1} -R {2} {3} -o {4}'
            ).format(
                "/bin/python3.5", self.quast_exe,
                self.ref_fasta,
                os.path.join(self.test_dir, "final_de_fere_novo_assembly",
                             "contigs.fasta"),
                os.path.join(self.test_dir, "quast_de_fere_novo")),
            str(
                '{0} {1} -R {2} {3} -o {4}'
            ).format(
                "/bin/python3.5", self.quast_exe,
                self.ref_fasta,
                os.path.join(self.test_dir, "final_de_novo_assembly",
                             "contigs.fasta"),
                os.path.join(self.test_dir, "quast_de_novo"))]
        for i, cmd in enumerate([x[0] for x in final_cmds]):
            self.assertEqual(final_spades_cmds_ref[i], cmd)
        # for i, cmd in enumerate([x[1] for x in final_cmds]):
        #     self.assertEqual(final_quast_cmds_ref[i], cmd)

    def test_def_decide_proceed_to_target_fail1(self):
        with self.assertRaises(ValueError):
            decide_proceed_to_target(target_len=25, logger=logger)
        with self.assertRaises(ValueError):
            decide_proceed_to_target(target_len=-1, logger=logger)
        with self.assertRaises(ValueError):
            decide_proceed_to_target(target_len=5.5, logger=logger)

    def test_def_decide_proceed_to_target(self):
        proceed = decide_proceed_to_target(target_len=4.5, logger=logger)
        dont_proceed = decide_proceed_to_target(target_len=None, logger=logger)
        self.assertTrue(proceed)
        self.assertFalse(dont_proceed)

    def test_samtools_depth_cmds_no_prep(self):
        """ check cmds: with region, no prep
        """
        self.assertEqual(
            (["samtools.iam sort test.bam > test_sorted.bam",
              "samtools.iam index test_sorted.bam"],
             "samtools.iam depth -r chrom1:33-99 test.bam"),
            make_samtools_depth_cmds(
                exe="samtools.iam",
                bam="test.bam", chrom="", start="", end="",
                region="chrom1:33-99", prep=False)
        )

    def test_samtools_depth_cmds_prep(self):
        """ check cmds: with region and prep
        """
        self.assertEqual(
            (["samtools.iam sort test.bam > test_sorted.bam",
              "samtools.iam index test_sorted.bam"],
             "samtools.iam depth -r chrom1:33-99 test_sorted.bam"),
            make_samtools_depth_cmds(
                exe="samtools.iam",
                bam="test.bam", chrom="", start="", end="",
                region="chrom1:33-99", prep=True)
        )

    def test_samtools_depth_cmds_coords(self):
        """ check cmds: without region, no prep
        """
        self.assertEqual(
            (["samtools.iam sort test.bam > test_sorted.bam",
              "samtools.iam index test_sorted.bam"],
             "samtools.iam depth -r chrom1:33-99 test.bam"),
            make_samtools_depth_cmds(
                exe="samtools.iam",
                bam="test.bam", chrom="chrom1", start=33, end=99,
                region=None, prep=False)
        )

    def test_make_modest_spades_cmd(self):
        """ test serialized allocation"""
        cmd = "spades.py --careful some args and stuff"
        self.assertEqual(
            "spades.py -t 4 -m 8 --careful some args and stuff",
            make_modest_spades_cmd(cmd=cmd, cores=4, memory=8,
                                   serialize=True, logger=logger))

    def test_make_modest_spades_cmd_parallel(self):
        """ test parallel allocation"""
        cmd = "spades.py --careful some args and stuff"
        self.assertEqual(
            "spades.py -t 1 -m 2 --careful some args and stuff",
            make_modest_spades_cmd(cmd=cmd, cores=4, memory=8,
                                   serialize=False, logger=logger))

    def test_make_modest_spades_cmd_split(self):
        """ test parallel allocation"""
        cmd = "spades.py --careful some args and stuff"
        self.assertEqual(
            "spades.py -t 2 -m 4 --careful some args and stuff",
            make_modest_spades_cmd(cmd=cmd, cores=4, memory=8, split=2,
                                   serialize=False, logger=logger))

    def test_make_modest_spades_cmd_odd(self):
        """ test parallel allocation, odd values"""
        cmd = "spades.py --careful some args and stuff"
        self.assertEqual(
            "spades.py -t 1 -m 3 --careful some args and stuff",
            make_modest_spades_cmd(cmd=cmd, cores=3, memory=11,
                                   serialize=False, logger=logger))

    def test_define_score_minimum_smalt(self):
        testargs = Namespace(mapper="smalt",
                             score_min=None)
        score_minimum = define_score_minimum(
            args=testargs, iteration=0,
            readlen=100, logger=logger)
        self.assertEqual(score_minimum, 50)

    def test_define_score_minimum_bwa(self):
        testargs = Namespace(mapper="bwa",
                             score_min=None)
        score_minimum = define_score_minimum(
            args=testargs, iteration=0,
            readlen=100, logger=logger)
        self.assertEqual(score_minimum, None)

    def test_define_score_minimum_explicit(self):
        testargs = Namespace(mapper="smalt",
                             score_min=33)
        score_minimum = define_score_minimum(
            args=testargs, iteration=0,
            readlen=100, logger=logger)
        self.assertEqual(score_minimum, 33)

    def test_define_score_minimum_badmapper(self):
        testargs = Namespace(mapper="cartography",
                             score_min=None)
        with self.assertRaises(AssertionError):
            score_minimum = define_score_minimum(
                args=testargs, iteration=0,
                readlen=100, logger=logger)

    def test_define_score_minimum_badmin(self):
        testargs = Namespace(mapper="bwa",
                             score_min=101)
        with self.assertRaises(ValueError):
            score_minimum = define_score_minimum(
                args=testargs, iteration=0,
                readlen=100, logger=logger)

    def test_bool_run_quast_false_wrongpython(self):
        with patch.object(sys, 'version_info') as v_info:
            v_info.minor = 6
            self.assertFalse(
                bool_run_quast("quast", logger))

    def test_bool_run_quast_false_none(self):
        self.assertFalse(bool_run_quast(None, logger))

    def test_bool_run_quast_true_regex(self):
        """ this just check the regex
        """
        pattern=r".* v(?P<version>[^\n,]+)"
        import re
        line = 1
        result = str.encode("QUAST v4.4\n")
        printout = result.decode("utf-8").split("\n")
        logger.warning(printout)
        this_version = None
        try:
            m = re.search(pattern, printout[line - 1])
        except IndexError as e:
            raise e
        print (m)
        if m:
            this_version = m.group('version').strip()
        logger.warning("this_version: %s", this_version)
        self.assertEqual("4.4", this_version)

    @unittest.skipIf(shutil.which("quast.py") is None and \
                     "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "Skipping this test on Travis CI. Too hard to debug")
    def test_bool_run_quast_true(self):
        self.assertTrue(bool_run_quast(shutil.which("quast"), logger))

    def test_make_quast_command(self):
        test_exes = Exes(python="/bin/python3.5",
                         samtools=self.samtools_exe,
                         quast=self.quast_exe,
                         smalt=self.smalt_exe,
                         spades=self.spades_exe,
                         bwa=self.bwa_exe,
                         check=False,
                         method="bwa")

        ref_cmd = str(
            "/bin/python3.5 quast ref.fasta /assembly/contigs.fasta -o " +
            "/here/quast_de_novo")
        test_cmd = make_quast_command(
            exes=test_exes, output_root="/here", ref="ref.fasta",
            assembly_subdir="/assembly/", name="de_novo", logger=logger)
        self.assertEqual(ref_cmd, test_cmd)

    def test_check_genbank_for_fasta(self):
        with self.assertRaises(ValueError):
            check_genbank_for_fasta(gb=self.ref_fasta, logger=logger)

    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "Skipping this test on Travis CI.")
    def test_filter_bam_AS(self):
        sourcebam = self.ref_bam_prefix + "_mapped.bam"
        filter_bam_AS(inbam=sourcebam,
                      outsam=os.path.join(self.test_dir, "filtered.sam"),
                      score=130, logger=logger)
        self.assertEqual(
            md5(os.path.join(self.test_dir, "filtered.sam")),
            md5(os.path.join(self.ref_dir, "test_bam_filtered_as130.sam")))
        self.to_be_removed.append(os.path.join(self.test_dir, "filtered.sam"))

    # @unittest.skipIf(shutil.which("samtools") is None or \
    #                  shutil.which("bcftools") is None,
    #                  "samtools executable not found, skipping." +
    #                  "If this isnt an error from travis deployment, you " +
    #                  "probably should install it")
    # def test_get_fasta_consensus_from_BAM(self):
    #     get_fasta_consensus_from_BAM(
    #         samtools_exe=shutil.which("samtools"),
    #         bcftools_exe=shutil.which("bcftools"),
    #         outfasta=os.path.join(self.test_dir, "bam_consensus.fasta"),
    #         ref=os.path.join(self.ref_dir, 'cluster1.fasta'),
    #         bam=self.ref_bam_prefix + "_mapped.bam",
    #         logger=logger)
    #     self.assertEqual(
    #         md5(os.path.join(self.test_dir, "bam_consensus.fasta")),
    #         md5(os.path.splitext(
    #             self.ref_fastq_with_iupac_bases)[0] + ".fasta"))

    def test_make_get_consensus_cmds_nofastq(self):
        cmd_list, outfasta = make_get_consensus_cmds(
            samtools_exe="samtools",
            bcftools_exe="bcftools",
            bam=self.ref_bam_prefix + "_mapped.bam",
            outfasta=None,
            ref=os.path.join(self.ref_dir, 'cluster1.fasta'),
            logger=logger)
        self.assertEqual(
            outfasta,
            self.ref_bam_prefix + "_mapped_consensus.fa")

    def test_make_get_consensus_cmds_withfastq(self):
        destfasta= self.ref_dir + "mapped.fasta"
        cmd_list, outfasta = make_get_consensus_cmds(
            samtools_exe="samtools",
            bcftools_exe="bcftools",
            bam=self.ref_bam_prefix + "_mapped.bam",
            outfasta=destfasta,
            ref=os.path.join(self.ref_dir, 'cluster1.fasta'),
            logger=logger)
        self.assertEqual(outfasta, destfasta)

    def test_make_get_consensus_cmds_check_cmds(self):
        ref_list = [
            "samtools faidx ref.fasta",
            "samtools sort bambam.bam > bambam_sorted.bam",
            "samtools index bambam_sorted.bam",
            "samtools mpileup -d8000 -EA -uf ref.fasta  bambam_sorted.bam | " +
            "bcftools call --ploidy 1 -c -Oz -o consensus.vcf.gz -",
            "tabix consensus.vcf.gz",
            "samtools faidx ref.fasta  | bcftools consensus " +
            "consensus.vcf.gz > out.fa"
        ]
        cmd_list, outfasta = make_get_consensus_cmds(
            samtools_exe="samtools",
            bcftools_exe="bcftools",
            bam="bambam.bam",
            outfasta="out.fa",
            old_method=True,
            ref="ref.fasta",
            logger=logger)
        for i in range(0, len(cmd_list)):
            self.assertEqual(cmd_list[i], ref_list[i])

    # def test_convert_fastq_to_fasta(self):
    #     outfasta = os.path.join(self.test_dir, "pass_fastq_to_fasta.fasta")
    #     convert_fastq_to_fasta(fastq=self.ref_fastq_with_iupac_bases,
    #                            outfasta=outfasta,
    #                            only_ATCG=True, logger=logger)
    #     self.assertEqual(
    #         md5(os.path.splitext(
    #             self.ref_fastq_with_iupac_bases)[0] + ".fasta"),
    #         md5(outfasta))

    def test_convert_fastq_to_fasta_failmulti(self):
        with self.assertRaises(ValueError):
            convert_fastq_to_fasta(
                fastq=self.ref_Rfastq,
                outfasta = os.path.join(
                    self.test_dir, "fail_fastq_to_fasta.fasta"),
                only_ATCG=True,
                logger=logger)

    def test_make_bwa_map_cmds_pe(self):
        ngs_ob = NgsLib(
            name="test",
            master=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)
        testmapping = LociMapping(
            name="test",
            iteration=1,
            assembly_subdir=self.test_dir,
            ref_fasta=self.ref_fasta,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        test_command_list = make_bwa_map_cmds(
            mapping_ob=testmapping, ngsLib=ngs_ob, cores=4,
            samtools_exe="samtools", bwa_exe="bwa",
            genome_fasta=self.ref_fasta,
            add_args='-L 0,0 -U 0 -a', logger=logger)
        ref_cmds = [
            "bwa index " + self.ref_fasta,  # index
            "bwa mem -t 4 -L 0,0 -U 0 -a -k 15 {0} {1} {2} | samtools view -bh - | samtools sort -o {3} - ".format(
                self.ref_fasta,  # 0
                self.ref_Ffastq,  # 1
                self.ref_Rfastq,  # 2
                testmapping.pe_map_bam),  # map
            "samtools view -bh {0} > {1}".format(
                testmapping.pe_map_bam,
                testmapping.mapped_bam_unfiltered)  # copy/merge
        ]

        for i, cmd in enumerate(test_command_list):
            self.assertEqual(cmd, ref_cmds[i])

    def test_make_bwa_map_cmds_s(self):
        ngs_ob = NgsLib(
            name="test",
            master=False,
            readS0=self.ref_Ffastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)
        testmapping = LociMapping(
            name="test",
            iteration=1,
            assembly_subdir=self.test_dir,
            ref_fasta=self.ref_fasta,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        test_command_list = make_bwa_map_cmds(
            mapping_ob=testmapping, ngsLib=ngs_ob, cores=4,
            samtools_exe="samtools", bwa_exe="bwa",
            genome_fasta=self.ref_fasta,
            add_args='-L 0,0 -U 0 -a', logger=logger)
        ref_cmds = [
            "bwa index " + self.ref_fasta,  # index
            "bwa mem -t 4 -L 0,0 -U 0 -a -k 15 {0} {1} | samtools view -bh - | samtools sort -o {2} - ".format(
                self.ref_fasta,  # 0
                self.ref_Ffastq,  # 1
                testmapping.s_map_bam),  # map
            "samtools view -bh {0} > {1}".format(
                testmapping.s_map_bam,
                testmapping.mapped_bam_unfiltered)  # copy/merge
        ]

        for i, cmd in enumerate(test_command_list):
            self.assertEqual(cmd, ref_cmds[i])

    def test_make_bwa_map_cmds_pes(self):
        ngs_ob = NgsLib(
            name="test",
            master=False,
            readS0=self.ref_Sfastq,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)
        testmapping = LociMapping(
            name="test",
            iteration=1,
            assembly_subdir=self.test_dir,
            ref_fasta=self.ref_fasta,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        test_command_list = make_bwa_map_cmds(
            mapping_ob=testmapping, ngsLib=ngs_ob, cores=4,
            samtools_exe="samtools", bwa_exe="bwa",
            genome_fasta=self.ref_fasta,
            add_args='-L 0,0 -U 0 -a', logger=logger)
        ref_cmds = [
            "bwa index " + self.ref_fasta,  # index
            "bwa mem -t 4 -L 0,0 -U 0 -a -k 15 {0} {1} {2} | samtools view -bh - | samtools sort -o {3} - ".format(
                self.ref_fasta,  # 0
                self.ref_Ffastq,  # 1
                self.ref_Rfastq,  # 2
                testmapping.pe_map_bam),  # map
            "bwa mem -t 4 -L 0,0 -U 0 -a -k 15 {0} {1} | samtools view -bh - | samtools sort -o {2} - ".format(
                self.ref_fasta,  # 0
                self.ref_Sfastq,  # 1
                testmapping.s_map_bam),  # map
            "samtools merge -f {0} {1} {2}".format(
                testmapping.mapped_bam_unfiltered,
                testmapping.pe_map_bam,
                testmapping.s_map_bam)  # copy/merge
        ]

        for i, cmd in enumerate(test_command_list):
            self.assertEqual(cmd, ref_cmds[i])

    def test_set_ref_as_contig_trust(self):
        self.assertEqual(
            set_ref_as_contig(ref_arg="trusted", map_percentage=90,
                              logger=logger),
            "trusted")
        self.assertEqual(
            set_ref_as_contig(ref_arg="infer", map_percentage=90,
                              logger=logger),

            "trusted")

    def test_set_ref_as_contig_untrust(self):
        self.assertEqual(
            set_ref_as_contig(ref_arg="untrusted", map_percentage=90,
                              logger=logger),
            "untrusted")
        self.assertEqual(
            set_ref_as_contig(ref_arg="infer", map_percentage=70,
                              logger=logger),
            "untrusted")

    def test_set_ref_as_contig_ignore(self):
        self.assertEqual(
            set_ref_as_contig(ref_arg="ignore", map_percentage=90,
                              logger=logger),
            None)

    def test_exclude_subassembly_based_on_coverage(self):
        gen = SeedGenome(
            max_iterations=1,
            genbank_path=self.ref_gb,
            clustered_loci_txt=self.test_loci_file,
            output_root=self.test_dir,
            logger=logger)
        # testmapping = LociMapping(
        #     name="test",
        #     iteration=1,
        #     assembly_subdir=self.test_dir,
        #     ref_fasta=self.ref_fasta,
        #     mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        gen.loci_clusters = parse_clustered_loci_file(
            filepath=gen.clustered_loci_txt,
            gb_filepath=gen.genbank_path,
            output_root=self.test_dir,
            padding=1000,
            circular=False,
            logger=logger)
        gen.loci_clusters[0].coverage_exclusion = True
        # self.assertEqual(
        #     exclude_subassembly_based_on_coverage(
        #         clu=gen.loci_clusters[0], iteration=1, logger=logger),
        #     0)
        self.assertEqual(
            exclude_subassembly_based_on_coverage(
                clu=gen.loci_clusters[0], iteration=0, logger=logger),
            2)
        # this should only be none or true, never false
        with self.assertRaises(AssertionError):
            gen.loci_clusters[0].coverage_exclusion = False
            exclude_subassembly_based_on_coverage(
                clu=gen.loci_clusters[0], iteration=0, logger=logger)
        gen.loci_clusters[0].coverage_exclusion = None
        self.assertEqual(None,
        exclude_subassembly_based_on_coverage(
                clu=gen.loci_clusters[0], iteration=0, logger=logger))

    def test_check_spades_extra_library_input(self):
        oddnum = "sg -g rm testdir everythin"
        # badcmd4 = "--mp3-1 path/to/file1 path/to/file1 --mp2-1 path/to/elsewhere"
        badcmd1 = "-f file1; rm testdir"
        badcmd2 = "sg -g s; || rm testdir"
        badcmd3 = "--pe3-1 test; rm testdir"
        badcmd4 = "--mp3-1 path/to/file1 -nsanopore path/to/file1 --mp2-1 path/to/elsewhere"
        badcmd5 = "--mp3-1 path/to/;cowsay;file1 -nsanopore path/to/file1 --mp2-1 path/to/elsewhere"

        with self.assertRaises(IndexError):
            check_spades_extra_library_input(oddnum)
        for badcmd in [badcmd1, badcmd2, badcmd3, badcmd4,badcmd5]:
            with self.assertRaises(ValueError):
                check_spades_extra_library_input(badcmd)
        goodcmd1 = "--pe2-1 path/to/files --pe2-1 path/to/elsewhere"
        goodcmd2 = "--pe9-1 path/to/files --nanopore path/to/elsewhere"
        goodcmd3 = "--mp3-1 path/to/files --mp2-1 path/to/elsewhere"
        goodcmd3 = "--s4 path/to/files --sanger path/to/elsewhere"
        for goodcmd in [goodcmd1, goodcmd2, goodcmd3]:
            check_spades_extra_library_input(goodcmd)


    def tearDown(self):
        """ delete temp files if no errors
        """

        for filename in self.to_be_removed:
            os.unlink(filename)
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
@unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                 "Skipping this test on Travis CI.")
class riboSeedDeep(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_riboSeed_tests")
        self.spades_dir = os.path.join(os.path.dirname(__file__),
                                       "output_riboSeed_tests",
                                       "SPAdes_results")
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.riboSeed_ref_dir = os.path.join(
            os.path.dirname(__file__), "references", "riboSeed_references")
        self.ref_gb = os.path.join(self.ref_dir,
                                   'NC_011751.1.gb')
        self.ref_fasta = os.path.join(self.test_dir,
                                      'cluster1.fasta')
        self.good_contig = os.path.join(self.ref_dir,
                                        'contigs.fasta')
        self.short_contig = os.path.join(self.ref_dir,
                                         'contigs.fasta')
        self.ref_Ffastq = os.path.join(self.ref_dir,
                                       'toy_reads1.fq')
        self.ref_Rfastq = os.path.join(self.ref_dir,
                                       'toy_reads2.fq')
        self.ref_Sfastq = os.path.join(self.ref_dir,
                                       'toy_readsS.fq')
        self.ref_bam_prefix = os.path.join(self.ref_dir,
                                           'test_bam_to_fastq')
        self.depthdir = os.path.join(self.ref_dir,
                                     "samtools_depth_test_files")
        self.sam = os.path.join(self.riboSeed_ref_dir,
                                "test_mapping.sam")
        self.samtxt = os.path.join(self.riboSeed_ref_dir,
                                   "some_sam_names.txt")

        self.smalt_exe = "smalt"
        self.bwa_exe = "bwa"
        self.samtools_exe = "samtools"
        self.spades_exe = "spades.py"
        self.quast_exe = "quast"
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
                              spades_exe="spades.py",
                              cores=2)
        self.cores = 2
        self.maxDiff = 2000
        self.to_be_removed = []
        self.startTime = time.time()
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)
        self.copy_fasta()

    def copy_fasta(self):
        """ make a disposable copy
        """
        shutil.copy(os.path.join(self.ref_dir, 'cluster1.fasta'),
                    self.ref_fasta)
        self.to_be_removed.append(self.ref_fasta)

    def test_map_with_bwa(self):
        # becasue multiple mapping are assingmed randomly (pseudorandomly?),
        # this accepts a range of expected results
        testmapping = LociMapping(
            name="test",
            iteration=1,
            assembly_subdir=self.test_dir,
            ref_fasta=self.ref_fasta,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        testngs = NgsLib(
            name="test",
            master=True,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.bwa_exe,
            logger=logger)

        map_to_genome_ref_bwa(
            mapping_ob=testmapping, ngsLib=testngs,
            genome_fasta=self.ref_fasta,
            cores=4, samtools_exe=self.samtools_exe,
            bwa_exe=self.bwa_exe, score_minimum=20,
            logger=logger)
        mapped_str = get_number_mapped(testmapping.pe_map_bam,
                                       samtools_exe=self.samtools_exe)
        nmapped = int(mapped_str[0:5])
        nperc = float(mapped_str[-13:-8])
        print("\nBWA: number of PE reads mapped: %f2" % nmapped)
        print("BWA: percentage of PE reads mapped: %2f" % nperc)
        ###
        testngs2 = NgsLib(
            name="test",
            master=True,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            readS0=self.ref_Sfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.bwa_exe,
            logger=logger)

        map_to_genome_ref_bwa(
            mapping_ob=testmapping, ngsLib=testngs2,
            genome_fasta=self.ref_fasta,
            cores=4, samtools_exe=self.samtools_exe,
            bwa_exe=self.bwa_exe, score_minimum=20,
            logger=logger)
        mapped_str2g = get_number_mapped(testmapping.pe_map_bam,
                                         samtools_exe=self.samtools_exe)
        mapped_str2 = get_number_mapped(testmapping.s_map_bam,
                                        samtools_exe=self.samtools_exe)
        nmapped2 = int(mapped_str2[0:5])
        nperc2 = float(mapped_str2[-13:-8])
        print("BWA: number of s reads mapped: %f2" % nmapped2)
        print("BWA: percentage of s reads mapped: %f2" % nperc2)

        ###
        self.assertTrue(13000 < nmapped < 14000)
        self.assertTrue(35.8 < nperc < 37.4)
        ###
        self.assertTrue(6300 < nmapped2 < 7000)
        self.assertTrue(35.0 < nperc2 < 36.8)

    def test_get_samtools_depths(self):
        """test bam generated by using SMALT to map the
        sample data reads (from the bam install check) to this reference
        (e coli genome with simplified chromosome name, gi12345)
        """
        print("starting depth test")
        test_bam = os.path.join(self.depthdir, "newref.bam")
        region = (1, 10000000)
        results = []
        for i in [region]:
            results.append(get_samtools_depths(
                samtools_exe=self.samtools_exe,
                bam=test_bam,
                chrom="gi12345", start=region[0],
                end=region[1],
                region=None,
                prep=True,
                logger=logger))
        self.assertEqual(round(results[0][1], 4), .9945)

    def test_fiddle_with_spades_exe(self):
        """ Verify the right python version is used for thisversion of SPAdes.
        This is a bad test, but there is no way that I know of to mock both
        subprocesses output and shutil.which and sys.executable.

        If it makes you feel better, I do have a hard time sleeping.
        """
        self.assertEqual(
            fiddle_with_spades_exe(shutil.which("spades.py"), logger=logger),
            sys.executable)

    def test_get_pysam_AS(self):
        """ test the getting of alignment scores from a sorted bam
        """
        test_bam_sorted = os.path.join(self.depthdir, "newref_sorted.bam")
        self.assertEqual(get_bam_AS(inbam=test_bam_sorted, logger=logger)[0:4],
                         [80, 57, 23, 89])

    def test_get_pysam_AS_unsorted(self):
        """ fail getting of alignment scores from an unsorted bam
        """
        test_bam = os.path.join(self.depthdir, "newref.bam")
        with self.assertRaises(ValueError):
            get_bam_AS(inbam=test_bam, logger=logger)

    def test_bam_to_fastq(self):
        testmapping = LociMapping(
            name="test",
            iteration=1,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        cmd, ngs = convert_bam_to_fastqs_cmd(mapping_ob=testmapping,
                                             ref_fasta="test_reference.fasta",
                                             samtools_exe=self.samtools_exe,
                                             which='mapped',
                                             single=True,
                                             source_ext="_bam",
                                             logger=logger)
        subprocess.run([cmd],
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
        self.assertTrue(os.path.getsize(ngs.readS0) > 0)


    def test_pysam_extract_reads(self):
        """ test extracting reads from a sam file
        """
        tempout = os.path.join(os.path.dirname(self.samtxt),
                               "temp_unmapped.sam")
        samref = os.path.join(os.path.dirname(self.samtxt),
                              "ref_unmapped.sam")
        pysam_extract_reads(
            sam=self.sam,
            textfile=self.samtxt,
            unmapped_sam=tempout, logger=logger)
        self.assertEqual(md5(tempout), md5(samref))
        self.to_be_removed.append(tempout)

    def test_prepare_next_mapping(self):
        gen = SeedGenome(
            max_iterations=1,
            genbank_path=self.ref_gb,
            clustered_loci_txt=self.test_loci_file,
            output_root=self.test_dir,
            logger=logger)
        testmapping = LociMapping(
            name="test",
            iteration=1,
            assembly_subdir=self.test_dir,
            ref_fasta=self.ref_fasta,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        testngs = NgsLib(
            name="test",
            master=True,
            make_dist=True,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe,
            logger=logger)
        gen.master_ngs_ob = testngs
        gen.loci_clusters = parse_clustered_loci_file(
            filepath=gen.clustered_loci_txt,
            gb_filepath=gen.genbank_path,
            output_root=self.test_dir,
            padding=1000,
            circular=False,
            logger=logger)
        add_coords_to_clusters(seedGenome=gen,
                               logger=logger)
        gen.next_reference_path = gen.ref_fasta
        #
        for cluster in gen.loci_clusters:
            cluster.master_ngs_ob = gen.master_ngs_ob

        map_to_genome_ref_smalt(
            mapping_ob=testmapping, ngsLib=testngs,
            genome_fasta=self.ref_fasta,
            cores=4, samtools_exe=self.samtools_exe,
            smalt_exe=self.smalt_exe, score_minimum=48,
            scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
            step=3, k=5, logger=logger)
        clu = gen.loci_clusters[0]
        self.assertEqual(len(clu.mappings), 0)
        prepare_next_mapping(
            cluster=clu, seedGenome=gen,
            flank=0,
            logger=logger)
        new_name = "NC_011751.1_cluster_{0}".format(clu.index)
        self.assertEqual(clu.mappings[-1].name, new_name)

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            os.unlink(filename)
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
@unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                 "Skipping this test on Travis CI.")
class riboSeedMain(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.args = Namespace(skip_contol=False, kmers="21,33,55,77,99",
                              spades_exe="spades.py",
                              cores=2)
        self.cores = 2
        self.maxDiff = 2000
        self.to_be_removed = []
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)
        self.copy_fasta()

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            os.unlink(filename)
        pass


if __name__ == '__main__':
    print("in a rush? run:  nosetests --with-coverage " +
          "tests/test_riboSeed.py:riboSeedShallow -v --cover-package")
    unittest.main()
