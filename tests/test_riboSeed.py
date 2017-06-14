# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas

"""
import sys
import logging
import shutil
# import subprocess
import os
import unittest
# import multiprocessing

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from argparse import Namespace

# I hate this line but it works :(
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "riboSeed"))


from pyutilsnrw.utils3_5 import md5, file_len, get_number_mapped

from riboSeed.riboSeed import SeedGenome, NgsLib, LociMapping, Exes, \
    map_to_genome_ref_smalt, map_to_genome_ref_bwa, \
    add_coords_to_clusters, partition_mapping, \
    convert_bam_to_fastqs_cmd, get_smalt_full_install_cmds,\
    generate_spades_cmd, estimate_distances_smalt, get_final_assemblies_cmds,\
    nonify_empty_lib_files, make_faux_genome, \
    evaluate_spades_success, prepare_next_mapping, make_mapped_partition_cmds,\
    make_unmapped_partition_cmds, make_quick_quast_table, \
    check_fastqs_len_equal, parse_subassembly_return_code, \
    make_spades_empty_check, get_samtools_depths, \
    decide_proceed_to_target, get_rec_from_generator, \
    check_kmer_vs_reads, make_samtools_depth_cmds, \
    parse_samtools_depth_results, make_modest_spades_cmd

from riboSeed.riboSnag import parse_clustered_loci_file, \
    extract_coords_from_locus, stitch_together_target_regions, \
    prepare_prank_cmd, prepare_mafft_cmd, \
    calc_Shannon_entropy, plot_scatter_with_anno, \
    profile_kmer_occurances, plot_pairwise_least_squares, make_msa


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
        self.ref_bam_prefix = os.path.join(self.ref_dir,
                                           'test_bam_to_fastq')
        self.smalt_exe = "smalt"
        self.bwa_exe = "bwa"
        self.samtools_exe = "samtools"
        self.spades_exe = "spades.py"
        self.quast_exe = "quast.py"
        self.python2_7_exe = "python2"
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
                              quast_exe="python2.7 quast.py",
                              cores=2)
        self.cores = 2
        self.maxDiff = 2000
        self.to_be_removed = []
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)
        self.copy_fasta()
        self.print_headsup()

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

    def test_Exes_bad_method(self):
        """check bad method arg given to Exes object"""
        with self.assertRaises(ValueError):
            Exes(samtools=self.samtools_exe,
                 quast=self.quast_exe,
                 smalt=self.smalt_exe,
                 python2_7=self.python2_7_exe,
                 spades=self.spades_exe,
                 bwa=self.bwa_exe,
                 method="bowtie")

    def test_Exes_bad_attribute(self):
        """ check bad instantiation of Exes object"""
        with self.assertRaises(AssertionError):
            Exes(samtools=None,
                 quast=self.quast_exe,
                 spades=self.spades_exe,
                 python2_7=self.python2_7_exe,
                 smalt=self.smalt_exe,
                 bwa=self.bwa_exe,
                 method="bwa")

    def test_Exes_bad_exe(self):
        """check Exes with nonexistant executable
        """
        with self.assertRaises(ValueError):
            Exes(samtools="nottheactualsamtools_exe",
                 quast=self.quast_exe,
                 python2_7=self.python2_7_exe,
                 smalt=self.smalt_exe,
                 spades=self.spades_exe,
                 bwa=self.bwa_exe,
                 method="bwa")

    def test_NgsLib(self):
        """ Can we create an NgsLib object correctly
        """
        # make a non-master object
        testlib_pe_s = NgsLib(
            name="test",
            master=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            readS0="dummy",
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)
        testlib_s = NgsLib(
            name="test",
            master=True,
            readF=None,
            readR=None,
            readS0=self.ref_Ffastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)
        self.assertEqual(testlib_s.libtype, "s_1")
        self.assertEqual(testlib_s.readlen, 145.0)
        self.assertEqual(testlib_s.liblist, [self.ref_Ffastq])
        # test unnamed fails
        with self.assertRaises(ValueError):
            NgsLib(
                name=None,
                master=False,
                readF=self.ref_Ffastq,
                readR=self.ref_Rfastq,
                readS0="dummy",
                ref_fasta=self.ref_fasta,
                mapper_exe=self.smalt_exe)
        self.assertEqual(testlib_pe_s.libtype, "pe_s")
        self.assertEqual(testlib_pe_s.readlen, None)
        # test fails with singe PE file
        with self.assertRaises(ValueError):
            NgsLib(
                name=None,
                master=False,
                readF=self.ref_Ffastq,
                readR=None,
                readS0="dummy",
                ref_fasta=self.ref_fasta,
                mapper_exe=self.smalt_exe)

        # check master files cannot bge deleted
        testlib_pe = NgsLib(
            name="test",
            master=True,
            make_dist=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe,
            logger=logger)
        self.assertEqual(
            1,  # return code for no deleting tool place
            testlib_pe.purge_old_files(master=testlib_pe_s, logger=logger))
        self.assertTrue(os.path.isfile(self.ref_Ffastq))
        self.assertTrue(os.path.isfile(self.ref_Rfastq))

        # test killer lib that tries to purge files that are in a master ob
        testlib_killer = NgsLib(
            name="test",
            master=False,
            make_dist=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe,
            logger=logger)
        self.assertEqual(
            1,  # return code for no deleting tool place
            testlib_killer.purge_old_files(master=testlib_pe_s, logger=logger))
        self.assertTrue(os.path.isfile(self.ref_Ffastq))
        self.assertTrue(os.path.isfile(self.ref_Rfastq))

    def test_SeedGenome(self):
        """ Can we create a seedGenome object
        """
        with open(self.ref_tiny_gb, "r") as gbf:
            rec = next(SeqIO.parse(gbf, "genbank"))
        gen = SeedGenome(
            max_iterations=2,
            clustered_loci_txt=self.test_loci_file,
            genbank_path=self.ref_tiny_gb,
            loci_clusters=None,
            output_root=self.test_dir)
        self.assertTrue(os.path.exists(os.path.join(self.test_dir,
                                                    "scannedScaffolds.fasta")))
        self.assertTrue(os.path.exists(self.test_dir))
        self.assertEqual(tuple(gen.seq_records)[0].id, rec.id)

    def test_SeedGenome_bad_instatiation(self):
        """ Does SeedGenome fail when missing attributes """
        with self.assertRaises(ValueError):
            SeedGenome(
                max_iterations=2,
                genbank_path=self.ref_gb,
                clustered_loci_txt=self.test_loci_file,
                loci_clusters=None,
                output_root=None)

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

    def test_add_coords_to_SeedGenome(self):
        """ can we parse a loci_coords file and add to seedGenome
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
        add_coords_to_clusters(seedGenome=gen, logger=logger)
        self.assertEqual(
            gen.loci_clusters[0].loci_list[0].start_coord, 4656045)
        self.assertEqual(
            gen.loci_clusters[0].loci_list[0].end_coord, 4657586)

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

    def test_generate_spades_cmds(self):
        """ can we make spades commands of various complexities
        """
        testmapping = LociMapping(
            name="test",
            iteration=1,
            assembly_subdir=self.test_dir,
            ref_fasta=self.ref_fasta,
            mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
        testngs1 = NgsLib(
            name="test",
            master=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)
        testngs2 = NgsLib(
            name="test",
            master=False,
            readF=self.ref_Ffastq,
            readR=self.ref_Rfastq,
            readS0=self.ref_Rfastq,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)
        # PE reads, prelim
        cmd1 = generate_spades_cmd(mapping_ob=testmapping, ngs_ob=testngs1,
                                   ref_as_contig='trusted',
                                   as_paired=True, addLibs="",
                                   prelim=True,
                                   k="21,33,55,77,99",
                                   spades_exe="spades.py",
                                   logger=logger)
        # PE with singletons
        cmd2 = generate_spades_cmd(mapping_ob=testmapping, ngs_ob=testngs2,
                                   ref_as_contig='trusted',
                                   as_paired=True, addLibs="",
                                   prelim=False,
                                   k="21,33,55,77,99",
                                   spades_exe="spades.py",
                                   logger=logger)
        # PE as sngle libraries
        cmd3 = generate_spades_cmd(mapping_ob=testmapping, ngs_ob=testngs1,
                                   ref_as_contig='trusted',
                                   as_paired=False, addLibs="",
                                   prelim=False,
                                   k="21,33,55,77,99",
                                   spades_exe="spades.py",
                                   logger=logger)
        # PE and singletons all as single libraries
        cmd4 = generate_spades_cmd(mapping_ob=testmapping, ngs_ob=testngs2,
                                   ref_as_contig='trusted',
                                   as_paired=False, addLibs="",
                                   prelim=False,
                                   k="21,33,55,77,99",
                                   spades_exe="spades.py",
                                   logger=logger)
        cmd1_ref = str(
            "spades.py --only-assembler --cov-cutoff off --sc --careful " +
            "-k 21,33,55,77,99 --pe1-1 {0} " +
            "--pe1-2 {1} --trusted-contigs {2}  -o {3}"
        ).format(
            self.ref_Ffastq, self.ref_Rfastq, self.ref_fasta, self.test_dir)
        self.assertEqual(cmd1, cmd1_ref)

        cmd2_ref = str(
            "spades.py --careful -k 21,33,55,77,99 --pe1-1 {0} " +
            "--pe1-2 {1} --pe1-s {2} --trusted-contigs {3}  -o {4}"
        ).format(
            self.ref_Ffastq, self.ref_Rfastq, self.ref_Rfastq,
            self.ref_fasta, self.test_dir)
        self.assertEqual(cmd2, cmd2_ref)

        cmd3_ref = str(
            "spades.py --careful -k 21,33,55,77,99 --pe1-s {0} " +
            "--pe2-s {1} --trusted-contigs {2}  -o {3}"
        ).format(
            self.ref_Ffastq, self.ref_Rfastq, self.ref_fasta, self.test_dir)
        self.assertEqual(cmd3, cmd3_ref)

        cmd4_ref = str(
            "spades.py --careful -k 21,33,55,77,99 --pe1-s {0} " +
            "--pe2-s {1} --pe3-s {1}  --trusted-contigs {2}  -o {3}"
        ).format(
            self.ref_Ffastq, self.ref_Rfastq, self.ref_fasta, self.test_dir)
        self.assertEqual(cmd4, cmd4_ref)

    def test_lib_check(self):
        """ does the NgsLib identify empty libraries
        """
        empty_file = os.path.join(self.test_dir, "test_not_real_file")
        # make an empty file
        with open(empty_file, 'w') as ef:
            pass
        ngs_ob = NgsLib(
            name="test",
            master=False,
            readF=self.ref_Ffastq,
            readR=empty_file,
            ref_fasta=self.ref_fasta,
            mapper_exe=self.smalt_exe)
        nonify_empty_lib_files(ngsLib=ngs_ob, logger=logger)
        self.assertTrue(ngs_ob.readR is None)
        self.to_be_removed.append(empty_file)

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
                                      final_contigs_dir=self.test_dir,
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
                                          final_contigs_dir=self.test_dir,
                                          skip_copy=True,
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
                                      final_contigs_dir=self.test_dir,
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
                                      final_contigs_dir=self.test_dir,
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
        """ Test all the possible combination of return codes givena particular
        spades results
        """
        shutil.copyfile(self.good_contig,
                        os.path.join(self.test_dir,
                                     "contigs.fasta"))
        # assert os.path.exists(self.ref_fasta), "WTF"
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
        add_coords_to_clusters(seedGenome=gen, logger=logger)
        clu = gen.loci_clusters[0]
        clu.mappings.append(testmapping)
        code_0 = evaluate_spades_success(
            clu,
            mapping_ob=testmapping,
            proceed_to_target=False, target_len=None,
            min_assembly_len=5000,
            read_len=300,
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
            read_len=300,
            min_assembly_len=10000,
            include_short_contigs=False, keep_best_contig=True,
            seqname='', logger=logger)
        self.assertEqual(code_2, 2)
        code_3 = evaluate_spades_success(
            clu,
            mapping_ob=testmapping2,
            proceed_to_target=False, target_len=None,
            read_len=300,
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
        # this mostly does system calls; cant really test smoothlu
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
        add_coords_to_clusters(seedGenome=gen, logger=logger)
        gen.loci_clusters[0].global_start_coord = 5000
        gen.loci_clusters[0].global_end_coord = 10000
        # print(gen.loci_clusters[0].__dict__)
        cmds, region = make_mapped_partition_cmds(
            cluster=gen.loci_clusters[0],
            mapping_ob=testmapping,
            seedGenome=gen,
            samtools_exe=self.samtools_exe,
            logger=logger)
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
        add_coords_to_clusters(seedGenome=gen, logger=logger)

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
        test_exes = Exes(samtools=self.samtools_exe,
                         quast=self.quast_exe,
                         python2_7=self.python2_7_exe,
                         smalt=self.smalt_exe,
                         spades=self.spades_exe,
                         bwa=self.bwa_exe,
                         check=False,
                         method="bwa")
        # I want this  generalized, so replace the acual path with spades.py
        test_exes.spades = "spades.py"
        test_exes.quast = "quast.py"
        test_exes.python2_7 = "python2"
        gen.ref_fasta = self.ref_fasta
        final_cmds, quast_reports = \
            get_final_assemblies_cmds(
                ref_as_contig="trusted",
                cores=4,
                serialize=True,
                memory=8,
                seedGenome=gen, exes=test_exes,
                skip_control=False, kmers="33,77,99", logger=logger)
        final_spades_cmds_ref = [
            str(
                "if [ -s {1} ] && [ -s {2} ] ; then " +
                "{0} -t 4 -m 8 --careful -k 33,77,99 --pe1-1 {1} " +
                "--pe1-2 {2} --trusted-contigs {3}  -o {4} " +
                "; else echo 'input lib not found, " +
                "skipping this SPAdes call' ; fi"
            ).format(
                self.spades_exe, self.ref_Ffastq, self.ref_Rfastq,
                gen.assembled_seeds,
                os.path.join(self.test_dir, "final_de_fere_novo_assembly")),
            str(
                "if [ -s {1} ] && [ -s {2} ] ; then " +
                "{0} -t 4 -m 8 --careful -k 33,77,99 --pe1-1 {1} " +
                "--pe1-2 {2}   -o {3} " +
                "; else echo 'input lib not found, " +
                "skipping this SPAdes call' ; fi"
            ).format(
                self.spades_exe, self.ref_Ffastq, self.ref_Rfastq,
                os.path.join(self.test_dir, "final_de_novo_assembly"))]
        final_quast_cmds_ref = [
            str(
                '{0} {1} -R {2} {3} -o {4}'
            ).format(
                self.python2_7_exe, self.quast_exe,
                self.ref_fasta,
                os.path.join(self.test_dir, "final_de_fere_novo_assembly",
                             "contigs.fasta"),
                os.path.join(self.test_dir, "quast_de_fere_novo")),
            str(
                '{0} {1} -R {2} {3} -o {4}'
            ).format(
                self.python2_7_exe, self.quast_exe,
                self.ref_fasta,
                os.path.join(self.test_dir, "final_de_novo_assembly",
                             "contigs.fasta"),
                os.path.join(self.test_dir, "quast_de_novo"))]
        for i, cmd in enumerate([x[0] for x in final_cmds]):
            self.assertEqual(final_spades_cmds_ref[i], cmd)
        for i, cmd in enumerate([x[1] for x in final_cmds]):
            self.assertEqual(final_quast_cmds_ref[i], cmd)

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

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            os.unlink(filename)
        pass


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class riboSeedDeep(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_riboseed_tests")
        self.spades_dir = os.path.join(os.path.dirname(__file__),
                                       "output_riboseed_tests",
                                       "SPAdes_results")
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
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
        self.ref_bam_prefix = os.path.join(self.ref_dir,
                                           'test_bam_to_fastq')
        self.smalt_exe = "smalt"
        self.bwa_exe = "bwa"
        self.samtools_exe = "samtools"
        self.spades_exe = "spades.py"
        self.quast_exe = "quast.py"
        self.python2_7_exe = "python2"
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
                              quast_exe="python2.7 quast.py",
                              cores=2)
        self.cores = 2
        self.maxDiff = 2000
        self.to_be_removed = []
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)
        self.copy_fasta()

    def copy_fasta(self):
        """ make a disposable copy
        """
        shutil.copy(os.path.join(self.ref_dir, 'cluster1.fasta'),
                    self.ref_fasta)
        self.to_be_removed.append(self.ref_fasta)

    @unittest.skipIf(shutil.which("bwa") is None or
                     shutil.which("quast.py") is None or
                     shutil.which("smalt") is None or
                     shutil.which("python2.7") is None or
                     shutil.which("spades.py") is None,
                     "bwa executable not found, skipping.If this isnt an " +
                     "error from travis deployment, you probably " +
                     "should install it")
    def test_Exes_(self):
        """check with  executable"""
        test_exes = Exes(samtools=self.samtools_exe,
                         quast=self.quast_exe,
                         python2_7=self.python2_7_exe,
                         smalt=self.smalt_exe,
                         spades=self.spades_exe,
                         bwa=self.bwa_exe,
                         method="bwa")
        self.assertEqual(test_exes.mapper, shutil.which(test_exes.bwa))

    # @unittest.skipIf(
    #     shutil.which("smalt") is None,
    #     "SMALT executable not found, skipping."+
    #     "If this isnt an error from travis deployment, you probably "+
    #     "should install it")
    # def test_estimate_distances_smalt(self):
    #     """ test estimate insert disances
    #     """
    #     if os.path.exists(self.test_estimation_file):
    #         print("warning! existing distance estimation file!")
    #     est_file = estimate_distances_smalt(outfile=self.test_estimation_file,
    #                                         smalt_exe=self.smalt_exe,
    #                                         ref_genome=self.ref_fasta,
    #                                         fastq1=self.ref_Ffastq,
    #                                         fastq2=self.ref_Rfastq,
    #                                         cores=1,
    #                                         logger=logger)
    #     ref_smi_md5 = "a444ccbcb486a8af29736028640e87cf"  # determined manually
    #     ref_sma_md5 = "4ce0c8b453f2bdabd73eaf8b5ee4f376"  # determined manually
    #     ref_mapping_len = 9271  # mapping doesnt have exact order, so cant md5
    #     self.assertEqual(ref_smi_md5, md5(str(est_file + ".smi")))
    #     self.assertEqual(ref_sma_md5, md5(str(est_file + ".sma")))
    #     self.assertEqual(ref_mapping_len, file_len(est_file))

    # @unittest.skipIf(
    #     shutil.which("smalt") is None,
    #     "smalt executable not found, skipping."+
    #     "If this isnt an error from travis deployment, you probably "+
    #     "should install it")
    # def test_map_with_smalt(self):
    #     # becasue multiple mapping are assingmed randomly (pseudorandomly?),
    #     # this accepts a range of expected results
    #     testmapping = LociMapping(
    #         name="test",
    #         iteration=1,
    #         assembly_subdir=self.test_dir,
    #         ref_fasta=self.ref_fasta,
    #         mapping_subdir=os.path.join(self.test_dir, "LociMapping"))
    #     testngs = NgsLib(
    #         name="test",
    #         master=True,
    #         make_dist=True,
    #         readF=self.ref_Ffastq,
    #         readR=self.ref_Rfastq,
    #         ref_fasta=self.ref_fasta,
    #         mapper_exe=self.smalt_exe,
    #         logger=logger)

    #     map_to_genome_ref_smalt(
    #         mapping_ob=testmapping, ngsLib=testngs,
    #         genome_fasta=self.ref_fasta,
    #         cores=4, samtools_exe=self.samtools_exe,
    #         smalt_exe=self.smalt_exe, score_minimum=48,
    #         scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
    #         step=3, k=5, logger=logger)
    #     mapped_str = get_number_mapped(testmapping.pe_map_bam,
    #                                    samtools_exe=self.samtools_exe)
    #     nmapped = int(mapped_str[0:5])
    #     nperc = float(mapped_str[-13:-8])
    #     print("\nSMALT: number of PE reads mapped: %f2" % nmapped)
    #     print("SMALT: percentage of PE reads mapped: %2f" % nperc)
    #     ###
    #     testngs2 = NgsLib(
    #         name="test",
    #         master=True,
    #         make_dist=True,
    #         readF=self.ref_Ffastq,
    #         readR=self.ref_Rfastq,
    #         readS0=self.ref_Rfastq,
    #         ref_fasta=self.ref_fasta,
    #         mapper_exe=self.smalt_exe,
    #         logger=logger)

    #     map_to_genome_ref_smalt(
    #         mapping_ob=testmapping, ngsLib=testngs2,
    #         genome_fasta=self.ref_fasta,
    #         cores=4, samtools_exe=self.samtools_exe,
    #         smalt_exe=self.smalt_exe, score_minimum=48,
    #         scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
    #         step=3, k=5, logger=logger)
    #     mapped_str2 = get_number_mapped(testmapping.pe_map_bam,
    #                                     samtools_exe=self.samtools_exe)
    #     mapped_str2 = get_number_mapped(testmapping.s_map_bam,
    #                                     samtools_exe=self.samtools_exe)
    #     nmapped2 = int(mapped_str2[0:5])
    #     nperc2 = float(mapped_str2[-13:-8])
    #     print("SMALT: number of s reads mapped: %f2" % nmapped2)
    #     print("SMALT: percentage of s reads mapped: %f2" % nperc2)

    #     ###
    #     self.assertTrue(12500 < nmapped < 12680)
    #     self.assertTrue(34.3 < nperc < 34.81)
    #     ###
    #     self.assertTrue(6400 < nmapped2 < 6500)
    #     self.assertTrue(34.3 < nperc2 < 35.8)

    @unittest.skipIf(shutil.which("bwa") is None,
                     "bwa executable not found, skipping." +
                     "If this isnt an error from travis deployment, you " +
                     "probably should install it")
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
            readS0=self.ref_Rfastq,
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

    @unittest.skipIf(shutil.which("samtools") is None,
                     "samtools executable not found, skipping." +
                     "If this isnt an error from travis deployment, you " +
                     "probably should install it")
    def test_get_samtools_depths(self):
        """test bam generated by using SMALT to map the
        sample data reads (from the bam install check) to this reference
        (e coli genome with simplified chromosome name, gi12345)
        """
        print("starting depth test")
        depthdir = os.path.join(self.ref_dir, "samtools_depth_test_files")
        # ref = os.path.join(depthdir, "test_ref.fasta")
        test_bam = os.path.join(depthdir, "newref.bam")
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
        # print(results)
        self.assertEqual(round(results[0][1], 4), .9945)

    @unittest.skipIf(shutil.which("smalt") is None, "smalt executable not found, skipping."+
                     "If this isnt an error from travis deployment, you probably "+
                     "should install it")
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
            cluster=clu, seedGenome=gen, samtools_exe=self.samtools_exe,
            flank=0,
            logger=logger)
        new_name = "NC_011751.1_cluster_{0}".format(clu.index)
        self.assertEqual(clu.mappings[-1].name, new_name)

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
