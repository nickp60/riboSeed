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
sys.dont_write_bytecode = True

from riboSeed.riboSnag import make_output_prefix, check_installed_tools,\
    copy_file, get_ave_read_len_from_fastq, get_number_mapped,\
    extract_mapped_and_mappedmates, keep_only_first_contig, md5,\
    combine_contigs, clean_temp_dir


def get_args():
    parser = argparse.ArgumentParser(
        description="test suite for pyutilsnrw repo")
    parser.add_argument("-k", "--keep_temps", dest='keep_temps',
                        action="store_true",
                        help="set if you want to inspect the output files",
                        default=False)
    args = parser.parse_args()
    return(args)


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class utils3_5TestCase(unittest.TestCase):
    # def init(self):
    #     pass
    # keep_temps = False
    # self.samtools = "test"
    # # samtools_exe = "samtoolz"
    # # def __init__(self, testname, keep_temps):
    # #     super(utils3_5TestCase, self).__init__(testname)
    # #     self.keep_temps = keep_temps
    # #     pass
    # print(self.samtools)
    def setUp(self):
        pass


def parse_clustered_loci_file(file):
    """
    """
    clusters={}
    index = 0
    try:
        with open(file, "r") as f:
            for line in f:
               clusters[index] = line.strip("\n").split(":")
               index = index + 1
    except:
        raise ValueError("problem parsing input clusters")
    return(clusters)


def extract_coords_from_locus(genome_seq_record, locus_tag_list=[],
                              verbose=True):
    """given a list of locus_tags, return a list of
    loc_number,coords, strand, product
    """
    loc_number = 0  # index for hits
    loc_list = []  # recipient structure
    for feat in genome_seq_record.features:
        try:
            if (feat.qualifiers.get("locus_tag")[0] in locus_tag_list) and\
               (feat.type == args.feature):
                #  SeqIO makes coords 0-based; the +1 below undoes that
                coords = [feat.location.start.position + 1,
                          feat.location.end.position]
                strand = feat.strand
                product = feat.qualifiers.get("product")
                locus_tag = feat.qualifiers.get("locus_tag")[0]
                loc_list.append([loc_number, coords, strand,
                                 product, locus_tag])
                loc_number = loc_number + 1
            else:
                pass
        except:
            pass
    if not loc_number > 0:
        raise ValueError("no hits found!")
        sys.exit(1)
    if verbose:
        print("Here are the detected region, strand, product, locus tag, \
               and subfeatures of the results:")
        pp.pprint(loc_list)
    return(loc_list)



def stitch_together_target_regions(genome_sequence, coords, flanking="500:500",
                                   within=50, minimum=50, replace=True,
                                   logger=None, verbose=True):
    """
    given a list from get_coords, usually of length 3 (16,5,and 23 rRNAs),
    return a string with the sequence of the region, replacing coding
    sequences with N's (or not, replace=False), and including the flanking
    regions upstream and down.

    revamped 20160913
    """
    if verbose and logger:
        log_status = logger.debug
    elif verbose:
        log_status = sys.stderr.write
    else:
        pass

    try:
        flank = [int(x) for x in flanking.split(":")]
        if len(flank) == 1:  # if only one value, use for both up and downstream
            flank.append(flank[0])
        assert(len(flank) == 2)
    except:
        raise ValueError("Error parsing flanking value; must either be " +
                         " integer or two colon-seapred integers")
    region = ''
    #TODO : make this safer
    smallest_feature = min([y[1] - y[0] for y in [ x[1] for x in coords]])
    # print(smallest_feature)
    if smallest_feature < (minimum):
        raise ValueError("invalid minimum! cannot exceed half of smallest " +
                         "feature, which is {0} in this case".format(
                             smallest_feature))
    # print("smallest featres")
    # print([y[1] - y[0] for y in [ x[1] for x in coords]])
    if verbose:
        for i in coords:
            log_status(str(i))
    #  This works as long as coords are never inreverse order
    global_start = min([y[0] for y in [x[1] for x in coords]]) - flank[0]
    global_end = max([y[1] for y in [x[1] for x in coords]]) + flank[1]
    full_seq = genome_sequence[global_start + 1: global_end]
    seq_with_ns = full_seq
    if replace:
        for i in coords:
            region_length = i[1][1] - i[1][0] - (2*within)
            # if dealing with short sequences
            if (i[1][1] - i[1][0]) < (2*within):
                # set within to retain minimum sequence length
                this_within = int((i[1][1] - i[1][0] - minimum) /2)
            else:
                # use default if not
                this_within = within
            rel_start = (i[1][0] + this_within) - global_start
            rel_end = (i[1][1] - this_within) - global_start
            # while rel_start < 0:
            #     region_within_start = int(region_within_start / 2)
            #     rel_start = (i[1][0] + region_within_start) - global_start
            # while rel_end < 0:
            #     region_within_end = int(region_within_end / 2)
            #     rel_end = (i[1][1] - region_within_end) - global_start
            seq_with_ns = str(seq_with_ns[0:rel_start] +
                              str("N"*region_length) +
                              seq_with_ns[rel_end: ])

    if verbose:
        log_status(str("exp length {0} \nact length {1}".format(
            global_end-global_start, len(full_seq))))
    if verbose:
        lb = 60
        for i in range(0, int(len(seq_with_ns) /lb)):
            log_status(str(full_seq[i*lb :lb+(i*lb)]+"\n"))
            log_status(str(seq_with_ns[i*lb :lb+(i*lb)]+"\n"))
            log_status("\n")
    return(seq_with_ns)















    def test_make_testing_dir(self):
        if not os.path.exists(testdirname):
            os.makedirs(testdirname)
        self.assertTrue(os.path.exists(testdirname))

    # def test_this_fails(self):
    #      self.assertEqual("pinecone", 42)

    def test_clean_temp_dir(self):
        """ I tried to do something like
        @unittest.skipUnless(clean_temp, "temporary files were retained")
        but couldnt get the variabel to be passed through.
        """
        if not os.path.exists(os.path.join(testdirname, "test_subdir")):
            os.makedirs(os.path.join(testdirname, "test_subdir"))
        clean_temp_dir(testdirname)

    def test_make_output_prefix(self):
        test_prefix = make_output_prefix(testdirname, "utils_3.5")
        self.assertEqual(test_prefix,
                         "".join([testdirname, os.path.sep, "utils_3.5"]))

    def test_check_installed_tools(self):
        """is pwd on all mac/linux systems?
        #TODO replace with better passing test
        """
        check_installed_tools(["pwd"])
        # test fails properly
        with self.assertRaises(SystemExit):
            check_installed_tools(["thisisnotapathtoanactualexecutable"])

    def test_md5_strings(self):
        """ minimal md5 examples
        """
        self.assertEqual(md5("thisstringisidenticalto", string=True),
                         md5("thisstringisidenticalto", string=True))
        self.assertNotEqual(md5("thisstringisntidenticalto", string=True),
                            md5("thisstringisnotidenticalto", string=True))

    def test_references_md5(self):
        """ is this paranoia, as well as bad testing?
        """
        test_pairs = [["3ba332f8a3b5d935ea6c4e410ccdf44b",
                       "references/combined_contigs_reference.fa"],
                      ["939fbf2c282091aec0dfa278b05e94ec",
                       "references/mapping_reference.bam"],
                      ["27944249bf064ba54576be83053e82b0",
                       "references/mapping_reference_mapped.sam"],
                      ["ac80c75f468011ba11e72ddee8560b33",
                       "references/md5_a.txt"],
                      ["ac80c75f468011ba11e72ddee8560b33",
                       "references/md5_b.txt"],
                      ["92fc8592819173343a75a40874d86144",
                       "references/md5_fail.txt"],
                      ["d6b0e5b28d0b4de431f10a03042ff37b",
                       "references/reads_reference.fastq"],
                      ["40ac496ec5b221636db81ce09e04c1d9",
                       "references/test_multiseqs_reference.fasta"],
                      ["920b5c9dc69fb2a9fed50b18f3e92895",
                       "references/test_only_first_reference.fasta"]]
        for i in test_pairs:
            self.assertEqual(i[0],
                             md5(os.path.join(os.path.dirname(__file__),
                                              i[1])))

    def test_md5_files(self):
        md5_a = md5(str(test_md5s_prefix + "_a.txt"))
        md5_b = md5(str(test_md5s_prefix + "_b.txt"))
        md5_fail = md5(str(test_md5s_prefix + "_fail.txt"))
        self.assertEqual(md5_a, md5_b)
        self.assertNotEqual(md5_a, md5_fail)

    def test_copy_file(self):
        if not os.path.exists(test_fastq_file):
            raise("test file is gone!  where is test_reads.fastq ?")
        new_path = copy_file(current_file=test_fastq_file,
                             dest_dir=testdirname,
                             name="newname.fastq", overwrite=False)
        # test path to copied file is constructed properly
        self.assertEqual(new_path, os.path.join(testdirname, "newname.fastq"))
        # test identity of files
        self.assertEqual(md5(new_path),
                         md5(os.path.join(testdirname, "newname.fastq")))
        # test overwrite exit
        with self.assertRaises(SystemExit):
            new_path = copy_file(current_file=test_fastq_file,
                                 dest_dir=testdirname,
                                 name="newname.fastq", overwrite=False)
        os.remove(new_path)

    def test_get_ave_read_len_from_fastq(self):
        """this probably could/should be refined to have a better test
        """
        mean_read_len = get_ave_read_len_from_fastq(test_fastq_file, N=5)
        self.assertEqual(217.8, mean_read_len)

    def test_get_number_mapped(self):
        """This is bad cause I hardcoded the path for samtools
        """
        result = get_number_mapped(test_bam_file, samtools_exe)
        reference = "151 + 0 mapped (0.56% : N/A)"
        self.assertEqual(result, reference)

    def test_extract_mapped_and_mappedmates(self):
        """ dont trust this if  make_output_prefix test fails
        some help from PSS on SO:
        http://stackoverflow.com/questions/16874598/
            how-do-i-calculate-the-md5-checksum-of-a-file-in-python
        """
        # copy files
        test_bam_dup = copy_file(current_file=test_bam_file,
                                 dest_dir=testdirname,
                                 name="", overwrite=False)
        test_mapped_dup = copy_file(current_file=test_sam_mapped_file,
                                    dest_dir=testdirname,
                                    name="", overwrite=False)
        ref_dir = os.path.join(testdirname)
        prefix = make_output_prefix(output_dir=ref_dir,
                                    name="mapping_reference")
        extract_mapped_and_mappedmates(map_results_prefix=prefix,
                                       fetch_mates=False,
                                       keep_unmapped=False,
                                       samtools_exe=samtools_exe)
        # reference mapping md5
        mapped_md5 = "27944249bf064ba54576be83053e82b0"
        md5_returned = md5(str(prefix + "_mapped.sam"))
        # Finally compare original MD5 with freshly calculated
        self.assertEqual(mapped_md5, md5_returned)
        # delete files created
        files_created = ["_mapped.bam",
                         "_mapped.bam.bai",
                         "_mapped.sam",
                         ".bam"]
        if mapped_md5 == md5_returned:
            for i in files_created:
                os.remove(str(prefix + i))

    def test_keep_only_first_contig(self):
        """copy_file
        """
        # copy to test dir
        copy_file(current_file=test_multifasta,
                  dest_dir=testdirname,
                  name='duplicated_multifasta.fasta', overwrite=False,
                  logger=None)
        path_to_dup = os.path.join(testdirname, "duplicated_multifasta.fasta")
        keep_only_first_contig(path_to_dup, newname="contig1")
        self.assertEqual(md5(path_to_dup), md5(test_singlefasta))
        # copy
        # copy_file(current_file=os.path.join(os.path.dirname(test_multifasta),
        #                                     "duplicated_multifasta.fasta"),
        #           dest_dir=os.path.dirname(test_multifasta),
        #           name=os.path.basename(test_multifasta), overwrite=True,
        #           logger=None)
        os.remove(path_to_dup)

    def test_combine_contigs(self):
        duplicated_multifasta = copy_file(current_file=test_multifasta,
                                          dest_dir=testdirname,
                                          name='multifasta_test_combine.fasta',
                                          overwrite=False,
                                          logger=None)
        for_first_contig = copy_file(current_file=test_multifasta,
                                     dest_dir=testdirname,
                                     name='single_fasta_test_combine.fasta',
                                     overwrite=False,
                                     logger=None)
        keep_only_first_contig(for_first_contig, newname="contig1")
        combined_contigs = combine_contigs(testdirname,
                                           pattern="*test_combine",
                                           contigs_name="combined_contigs.fa",
                                           ext=".fasta",
                                           verbose=True)
        self.assertEqual(md5(test_combined), md5(combined_contigs))
        for i in [duplicated_multifasta, for_first_contig, combined_contigs]:
            os.remove(i)

    def tearDown(self):
        pass
##  Functions to write tests for


# def setup_protein_blast(input_file, input_type="fasta", dbtype="prot",
#                         title="blastdb", out="blastdb",
#                         makeblastdb_exe=''):
#     """
#     This runs make blast db with the given parameters
#     requires logging, os, subprocess, shutil
#     """
#     logger = logging.getLogger(__name__)
#     #logging.getLogger(name=None)
#     logger.debug("TESTING I 2 3!")
#     if makeblastdb_exe == '':
#         makeblastdb_exe = shutil.which("makeblastdb")
#     makedbcmd = str("{0} -in {1} -input_type {2} -dbtype {3} " +
#                     "-title {4} -out {5}").format(makeblastdb_exe,
#                         input_file, input_type, dbtype, title, out)
#     logger.info("Making blast db: {0}".format(makedbcmd))
#     try:
#         subprocess.run(makedbcmd, shell=sys.platform != "win32",
#                        stdout=subprocess.PIPE,
#                        stderr=subprocess.PIPE, check=True)
#         logging.debug("BLAST database '{0}' created here: {1}".format(
#             title, out))
#         return(0)
#     except:
#         logging.error("Something bad happened when trying to make " +
#                       "a blast database")
#         sys.exit(1)


# def run_blastp(input_file, database_name, outfmt, blastp_exe=''):
#     """
#     requires logging subprocess, os, shutil
#     """
#     #logger = logging.getLogger(name=None)
#     logger = logging.getLogger(__name__)
#     output_file = os.path.join(os.path.split(input_file)[0],
#                                str(os.path.splitext(
#                                    os.path.basename(input_file))[0] +
#                                    "_blast_hits.tab"))
#     if blastp_exe == '':
#         blastp_exe = shutil.which("blastp")
#     blastpcmd = str("{0} -db {1} -query {2} -out {3} -outfmt " +
#                     "{4}").format(blastp_exe, database_name, input_file,
#                                   output_file, outfmt)
#     logger.info("Running blastp: {0}".format(blastpcmd))
#     try:
#         subprocess.run(blastpcmd, shell=sys.platform != "win32",
#                        stdout=subprocess.PIPE,
#                        stderr=subprocess.PIPE, check=True)
#         logger.debug("Results from BLASTing {0} are here: {1}".format(
#             input_file, output_file))
#         return(0)
#     except:
#         logger.error("Something bad happened when running blast")
#         sys.exit(1)


# #def merge_outfiles(filelist, outfile_name):
# def merge_blast_tab_outfiles(filelist, outfile_name):
#     """
#     #TODO needs a test for headers
#     for combining tab formated blast output format 6
#     returns 0 if successful
#     requires logging
#     """
#     # only grab .tab files, ie, the blast output
#     logger = logging.getLogger()
#     filelist = [i for i in filelist if i.split(".")[-1:] == ['tab']]
#     if len(filelist) == 1:
#         logger.warning("only one file found! no merging needed")
#         return(0)
#     elif len(filelist) == 0:
#         logger.error("filelist empt; cannot perform merge!")
#         return(1)
#     else:
#         logger.info("merging all the blast results to %s" % outfile_name)
#         nfiles = len(filelist)
#         fout = open(outfile_name, "a")
#         # first file:
#         for line in open(filelist[0]):
#             fout.write(line)
#         #  now the rest:
#         for num in range(1, nfiles):
#             f = open(filelist[num])
#             for line in f:
#                 fout.write(line)
#             f.close()  # not really needed
#         fout.close()
#         return(0)


# def cleanup_output_to_csv(infile,
#                           accession_pattern='(?P<accession>[A-Z _\d]*\.\d*)'):
#     """
#     given .tab from merge_blast_tab_outfiles, assign pretty column names,
#     """
#     logger = logging.getLogger(name=None)
#     print("cleaning up the csv output")
#     colnames = ["query_id", "subject_id", "identity_perc", "alignment_length", "mismatches",
#                 "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
#     csv_results = pd.read_csv(open(infile), comment="#", sep="\t", names=colnames)
#     #This regex will probably break things rather badly before too long...
#     # it looks for capital letter and numbers, dot, number, ie SHH11555JJ8.99
#     csv_results["accession"] = csv_results.query_id.str.extract(accession_pattern)
#     # write out results with new headers or with new headers and merged metadat from accessions.tab
#     genes = open(genelist, "r")
#     genedf = pd.read_csv(genes, sep=",")
#     output_path_csv = str(os.path.splitext(infile)[0]+".csv")
#     results_annotated = pd.merge(csv_results, genedf, how="left",  on="accession")
#     results_annotated.to_csv(open(output_path_csv, "w"))
#     print("wrote final csv to %s" % output_path_csv)
# #%%


if __name__ == '__main__':
    args = get_args()
    curdir = os.getcwd()
    # samtools_exe = args.samtools_exe
    testdirname = os.path.join(os.path.dirname(__file__),
                               "output_utils3_5_tests")
    test_fastq_file = os.path.join(os.path.dirname(__file__),
                                   str("references" + os.path.sep +
                                       'reads_reference.fastq'))
    test_bam_file = os.path.join(os.path.dirname(__file__),
                                 str("references" + os.path.sep +
                                     "mapping_reference.bam"))
    test_sam_mapped_file = os.path.join(os.path.dirname(__file__),
                                        str("references" + os.path.sep +
                                            "mapping_reference_mapped.sam"))
    test_multifasta = os.path.join(os.path.dirname(__file__),
                                   str("references" + os.path.sep +
                                       "test_multiseqs_reference.fasta"))
    test_singlefasta = os.path.join(os.path.dirname(__file__),
                                    str("references" + os.path.sep +
                                        "test_only_first_reference.fasta"))
    test_combined = os.path.join(os.path.dirname(__file__),
                                 str("references" + os.path.sep +
                                     "combined_contigs_reference.fa"))
    test_md5s_prefix = os.path.join(os.path.dirname(__file__),
                                    str("references" + os.path.sep +
                                        "md5"))
    # utils3_5TestCase.keep_temps = args.keep_temps
    samtools_exe = "samtools"
    unittest.main()
    # suite = unittest.TestSuite()
    # keep_temps=args.keep_temps
    # suite = unittest.defaultTestLoader.loadTestsFromTestCase(utils3_5TestCase)
    # suite.addTest(utils3_5TestCase("tearDown", keep_temps))
    # unittest.utils3_5TestCase().run(suite)
    # if not args.keep_temps:
    #     os.rmdir(testdirname)
