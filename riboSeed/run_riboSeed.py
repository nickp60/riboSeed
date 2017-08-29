#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import pkg_resources
import sys
import os
import shutil
import subprocess
def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="Given a directory of one or more chromosomes as fasta " +
        "files, this facilitates reannotation of rDNA regions with Barrnap " +
        " and outputs all sequences as a single, annotated genbank file",
        add_help=False)  # to allow for custom help
    parser.add_argument("contigs", action="store",
                        help="either a (multi)fasta or a directory " +
                        "containing one or more chromosomal " +
                        "sequences in fasta format")

    requiredNamed = parser.add_argument_group('required named arguments')
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-c", "--config", dest='config', action="store",
                          help="config file; if none given, create one; " +
                          "default: %(default)s", default=os.getcwd(),
                          type=str, required=True)
    optional.add_argument("-o", "--output", dest='output', action="store",
                          help="output directory; " +
                          "default: %(default)s",
                          default=os.path.join(
                              os.getcwd(),
                              str(time.strftime("%Y-%m-%dT%H:%M") +
                                  "_riboSeed_pipeline_results"), ""),
                          type=str)

    optional.add_argument("-k", "--kingdom", dest='kingdom',
                          action="store",
                          choices=["bac", "euk", "arc", "mito"],
                          help="whether to look for eukaryotic, archaeal, or" +
                          " bacterial rDNA; " +
                          "default: %(default)s", default="bac",
                          type=str)
    optional.add_argument("-t", "--id_thresh", dest='id_thresh',
                          action="store", type=float,
                          help="partial rRNA hits below this threshold will " +
                          "be ignored. " +
                          "default: %(default)s", default=0.5,
                          )
    optional.add_argument("-n", "--name", dest='name',
                          action="store",
                          # default="contig",
                          help="name to give the contig files; "
                          "default: infered from file",
                          type=str)
    optional.add_argument("-s", "--specific_features",
                          help="colon:separated -- specific features"
                          "; default: %(default)s",
                          default='16S:23S:5S', type=str)
    # TODO Implement this for work with really fragmented genomes
    # optional.add_argument("--nocluster",
    #                       help="do not bother clustering; treat all " +
    #                       "occurances on a sequence as a cluster" +\
    #                       "default: %(default)s", action='store_true',
    #                       default=False, dest="nocluster")
    optional.add_argument("-c", "--clusters",
                          help="number of rDNA clusters;"
                          "if submitting multiple records, must be a "
                          "colon:separated list whose length matches number "
                          "of genbank records.  Default is inferred from "
                          "specific feature with fewest hits", default='',
                          type=str, dest="clusters")
    optional.add_argument("-c", "--config", dest='config', action="store",
                          help="config file; if none given, use the default; " +
                          "default: %(default)s", default=os.getcwd(),
                          type=str, required=True)
    optional.add_argument("-F", "--fastq1", dest='fastq1', action="store",
                          help="forward fastq reads, can be compressed",
                          type=str, default=None)
    optional.add_argument("-R", "--fastq2", dest='fastq2', action="store",
                          help="reverse fastq reads, can be compressed",
                          type=str, default=None)
    optional.add_argument("-S1", "--fastq_single1", dest='fastqS1',
                          action="store",
                          help="single fastq reads", type=str, default=None)

    optional.add_argument("-j", "--just_seed", dest='just_seed',
                          action="store_true",
                          default=False,
                          help="Don't do an assembly, just generate the long" +
                          " read 'seeds'; default: %(default)s")
    optional.add_argument("-n", "--experiment_name", dest='exp_name',
                          action="store",
                          help="prefix for results files; " +
                          "default: %(default)s",
                          default="riboSeed", type=str)
    optional.add_argument("-l", "--flanking_length",
                          help="length of flanking regions, in bp; " +
                          "default: %(default)s",
                          default=1000, type=int, dest="flanking")
    optional.add_argument("-c", "--cores", dest='cores', action="store",
                          default=None, type=int,
                          help="cores used" +
                          "; default: %(default)s")
    optional.add_argument("--memory", dest='memory', action="store",
                          default=8, type=int,
                          help="cores for multiprocessing" +
                          "; default: %(default)s")
    optional.add_argument("-k", "--kmers", dest='kmers', action="store",
                          default="21,33,55,77,99,127", type=str,
                          help="kmers used for final assembly" +
                          ", separated by commas such as" +
                          "21,33,55,77,99,127 . Can be set to 'auto', where " +
                          "SPAdes chooses.  We ensure kmers are not " +
                          "too big or too close to read length" +
                          "; default: %(default)s")
    optional.add_argument("-p", "--pre_kmers", dest='pre_kmers',
                          action="store",
                          default="21,33,55,77,99", type=str,
                          help="kmers used during seeding assemblies, " +
                          "separated bt commas" +
                          "; default: %(default)s")
    optional.add_argument("-d", "--min_flank_depth",
                          help="a subassembly will not be performed if this " +
                          "minimum depth is not achieved on both the 3' and" +
                          "5' end of the pseudocontig. " +
                          "default: %(default)s",
                          default=0, dest="min_flank_depth", type=float)
    optional.add_argument("--clean_temps", dest='clean_temps',
                          default=False, action="store_true",
                          help="if --clean_temps, mapping files will be " +
                          "removed once they are no no longer needed during " +
                          "the mapping iterations to save space; " +
                          "default: %(default)s")
    optional.add_argument("-i", "--iterations", dest='iterations',
                          action="store",
                          default=3, type=int,
                          help="if iterations>1, multiple seedings will " +
                          "occur after subassembly of seed regions; " +
                          "if setting --target_len, seedings will continue " +
                          "until --iterations are completed or --target_len"
                          " is matched or exceeded; " +
                          "default: %(default)s")
    optional.add_argument("-v", "--verbosity", dest='verbosity',
                          action="store",
                          default=2, type=int, choices=[1, 2, 3, 4, 5],
                          help="Logger writes debug to file in output dir; " +
                          "this sets verbosity level sent to stderr. " +
                          " 1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
                          "default: %(default)s")
    optional.add_argument("-t", "--threads", dest='threads',
                          action="store",
                          default=1, type=int,
                          choices=[1, 2, 4],
                          help="if your cores are hyperthreaded, set number " +
                          "threads to the number of threads per processer." +
                          "If unsure, see 'cat /proc/cpuinfo' under 'cpu " +
                          "cores', or 'lscpu' under 'Thread(s) per core'." +
                          ": %(default)s")
    optional.add_argument("-z", "--serialize", dest='serialize',
                          action="store_true",
                          default=False,
                          help="if --serialize, runs seeding and assembly " +
                          "without multiprocessing. This is recommended for " +
                          "machines with less than 8GB RAM: %(default)s")
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")


def main():
    root_dir = os.path.join(os.getcwd(), "integration_test_results")
    os.makedirs(root_dir, exist_ok=True)
    cmds = [
        # riboScan
        "{0} {1} {2}{3} -o {4} -k bac".format(
            sys.executable,
            shutil.which("riboScan.py"),
            os.path.join(os.path.dirname(fasta), ""),
            os.path.basename(fasta), os.path.join(root_dir, "scan", "")),
        # riboSelect
        "{0} {1} {2}scannedScaffolds.gb -o {3}".format(
        sys.executable,
        shutil.which("riboSelect.py"),
        os.path.join(root_dir, "scan", ""),
        os.path.join(root_dir, "select","")),
        # riboSnag
        "{0} {1} {2}scannedScaffolds.gb {3}riboSelect_grouped_loci.txt -o {4} ".format(
            sys.executable,
            shutil.which("riboSnag.py"),
            os.path.join(root_dir, "scan", ""),
            os.path.join(root_dir, "select", ""),
            os.path.join(root_dir, "snag", "")),
        # riboSeed
        "{0} {1} -r {2}scannedScaffolds.gb {3}riboSelect_grouped_loci.txt -o {4} -F {5} -R {6} --serialize -v 1 ".format(
            sys.executable,
            shutil.which("riboSeed.py"),
            os.path.join(root_dir, "scan", ""),
            os.path.join(root_dir, "select", ""),
            os.path.join(root_dir, "seed", ""),
            fastq1,
            fastq2)
    ]
    for i in cmds:
        print(i)
        subprocess.run([i],
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    print("finished running integration test!")





    import importlib.util
    spec = importlib.util.spec_from_file_location("test", "./test.py")
    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    print(config.BARRNAP_EXE)



#  needed args
config_file
contigs (scan)
outdir (none)
kingdom (scan)
specific_features (select)
clusters (select)
reads (seed)
just_seed (seed)
name (seed)
flanking_length
cores
threads
memory
serialize
kmers
prekemers
min_flanking_depth
clean temps
iterationso
# config
id_trhesh (scan)
contig_name (scan)
min_length (scan)
feature (select)
debug (select)
map_method
score_min
min_assembly_length
include_shorts
subtract
clean temps
skip_control
ref_as_contig
target_len
mapper_args
smalt_scoring
