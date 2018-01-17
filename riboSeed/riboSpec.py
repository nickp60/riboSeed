#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# Copyright 2017, National University of Ireland and The James Hutton Insitute
# Author: Nicholas Waters
#
# This code is part of the riboSeed package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""
Created on Wed Jan 17 14:44:30 2018

We have developed a scheme by which to predict the number of rDNAs in a
 genome based on an assembly graph created by spades.
This module is designed to


"""

import sys
import os
import re

class FastgNode(object):
    def __init__(self, name=None, length=None, cov=None,
                 reverse_complimented=None, neighbor_list=None):
        self.name = name
        self.length = length
        self.cov = cov
        self.neighbor_list = neighbor_list
        self.reverse_complimented = reverse_complimented
    def __str__(self):
        return str("Node: {0}\nNeighbors: {1}\nLength: {2}\nCoverage: " +
                   "{3}\nReverse_Complimented?: {4}").format(
                       self.name,
                       "None" if self.neighbor_list is None else
                           ",".join([str(x.name) for x in self.neighbor_list]),
                       "None" if self.length is None else self.length,
                       "None" if self.cov is None else self.cov,
                       str(self.reverse_complimented)
                   )

class FastgNodeNeighbor(object):
    def __init__(self, name=None,reverse_complimented=None):
        self.name = name
        self.reverse_complimented = reverse_complimented
    def __str__(self):
        return "NeighborNode: {0}\nReverse_Complimented?: {1}".format(
            self.name,
            "Yes" if self.reverse_complimented else "No"
        )


def get_args():  # pragma: no cover
    """
    """
    parser = argparse.ArgumentParser(prog="ribo spec",
        description="Given either an assembly graph or a mapping file, " +
        "determine whether the number of rDNAs appears to match the " +
        "reference",
        add_help=False)  # to allow for custom help
    parser.add_argument("clustered_loci_txt", action="store",
                        help="output from riboSelect")
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-r", "--reference_genbank",
                               dest='reference_genbank',
                               action="store", default='', type=str,
                               help="genbank reference, used to estimate " +
                               "insert sizes, and compare with QUAST",
                               required=True)
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-F", "--fastq1", dest='fastq1', action="store",
                          help="forward fastq reads, can be compressed",
                          type=str, default=None)
    optional.add_argument("-R", "--fastq2", dest='fastq2', action="store",
                          help="reverse fastq reads, can be compressed",
                          type=str, default=None)
    optional.add_argument("-S1", "--fastq_single1", dest='fastqS1',
                          action="store",
                          help="single fastq reads", type=str, default=None)
    # parameters for run
    optional.add_argument("-l", "--flanking_length",
                          help="length of flanking regions, in bp; " +
                          "default: %(default)s",
                          default=1000, type=int, dest="flanking")
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
    optional.add_argument("--mapper", dest='method',
                          action="store", choices=["smalt", "bwa"],
                          help="available mappers: smalt and bwa; " +
                          "default: %(default)s",
                          default='bwa', type=str)
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
    optional.add_argument("--force_kmers", dest="force_kmers",
                          action="store_true",
                          default=False,
                          help="skip checking to see if kmerchoice is " +
                          "appropriate to read length. Sometimes kmers " +
                          "longer than reads can help in the final assembly," +
                          " as the long reads generated by riboSeed contain " +
                          "kmers longer than the read length")
    optional.add_argument("-s", "--score_min", dest='score_min',
                          action="store",
                          default=None, type=int,
                          help="If using smalt, this sets the '-m' param; " +
                          "default with smalt is inferred from " +
                          "read length. If using BWA, reads mapping with AS" +
                          "score lower than this will be rejected" +
                          "; default with BWA is half of read length")
    optional.add_argument("-a", "--min_assembly_len", dest='min_assembly_len',
                          action="store",
                          default=6000, type=int,
                          help="if initial SPAdes assembly largest contig " +
                          "is not at least as long as --min_assembly_len, " +
                          "reject. Set this to the length of the seed " +
                          "sequence; if it is not achieved, seeding across " +
                          "regions will likely fail; default: %(default)s")
    optional.add_argument("--include_shorts", dest='include_short_contigs',
                          action="store_true",
                          default=False,
                          help="if assembled contig is smaller than  " +
                          "--min_assembly_len, contig will still be included" +
                          " in assembly; default: inferred")
    optional.add_argument("--damn_the_torpedos", dest='damn_the_torpedos',
                          action="store_true",
                          default=False,
                          help="Ignore certain errors, full speed ahead!")
    optional.add_argument("--subtract", dest='subtract',
                          action="store_true",
                          default=False,
                          help="if --subtract reads already used in previous" +
                          "round of subassembly will not be included in " +
                          "subsequent rounds.  This can lead to problems " +
                          "with multiple mapping and inflated coverage.")
    optional.add_argument("--linear",
                          help="if genome is known to not be circular and " +
                          "a region of interest (including flanking bits) " +
                          "extends past chromosome end, this extends the " +
                          "seqence past chromosome origin forward by " +
                          "--padding; " +
                          "default: %(default)s",
                          default=False, dest="linear", action="store_true")
    optional.add_argument("-d", "--min_flank_depth",
                          help="a subassembly will not be performed if this " +
                          "minimum depth is not achieved on both the 3' and" +
                          "5' end of the pseudocontig. " +
                          "default: %(default)s",
                          default=0, dest="min_flank_depth", type=float)
    optional.add_argument("--ref_as_contig", dest='ref_as_contig',
                          action="store", type=str,
                          default="infer",
                          choices=["ignore", "infer", "trusted", "untrusted"],
                          help="ignore: reference will not be used in " +
                          "subassembly. trusted: SPAdes will use the seed" +
                          " sequences as a --trusted-contig; untrusted: " +
                          "SPAdes will treat as --untrusted-contig. " +
                          "infer: if mapping percentage " +
                          "over 80%%, 'trusted'; else 'untrusted'." +
                          " See SPAdes docs for details.  default: infer")
    optional.add_argument("--clean_temps", dest='clean_temps',
                          default=False, action="store_true",
                          help="if --clean_temps, mapping files will be " +
                          "removed once they are no no longer needed during " +
                          "the mapping iterations to save space; " +
                          "default: %(default)s")
    optional.add_argument("--skip_control", dest='skip_control',
                          action="store_true",
                          default=False,
                          help="if --skip_control, no de novo " +
                          "assembly will be done; default: %(default)s")
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
    optional.add_argument("--target_len", dest='target_len', action="store",
                          default=None, type=float,
                          help="if set, iterations will continue until " +
                          "contigs reach this length, or max iterations (" +
                          "set by --iterations) have been completed. Set as " +
                          "fraction of original seed length by giving a " +
                          "decimal between 0 and 5, or set as an absolute " +
                          "number of base pairs by giving an integer greater" +
                          " than 50. Not used by default")
    optional.add_argument("-z", "--serialize", dest='serialize',
                          action="store_true",
                          default=False,
                          help="if --serialize, runs seeding and assembly " +
                          "without multiprocessing. This is recommended for " +
                          "machines with less than 8GB RAM: %(default)s")
    optional.add_argument("--consensus", dest='initial_consensus',
                          action="store_true",
                          default=False,
                          help="if --initial_consensus, " +
                          "generate a mpileup-based consesnsus instead of " +
                          "doing a proper spades subassembly")
    optional.add_argument("--smalt_scoring", dest='smalt_scoring',
                          action="store",
                          default="match=1,subst=-4,gapopen=-4,gapext=-3",
                          help="if mapping with SMALT, " +
                          "submit custom smalt scoring via smalt -S " +
                          "scorespec option; default: %(default)s")
    optional.add_argument("--mapper_args", dest='mapper_args',
                          action="store",
                          default="-L 0,0 -U 0 -a",
                          help="submit custom parameters to mapper. " +
                          "And by mapper, I mean bwa, cause we dont support " +
                          "this option for SMALT, sorry. " +
                          "This requires knowledge of your chosen mapper's " +
                          "optional arguments. Proceed with caution!  " +
                          "default: %(default)s")
    # # TODO  Make these check a config file
    optional.add_argument("--spades_exe", dest="spades_exe",
                          action="store", default="spades.py",
                          help="Path to SPAdes executable; " +
                          "default: %(default)s")
    optional.add_argument("--samtools_exe", dest="samtools_exe",
                          action="store", default="samtools",
                          help="Path to samtools executable; " +
                          "default: %(default)s")
    optional.add_argument("--smalt_exe", dest="smalt_exe",
                          action="store", default="smalt",
                          help="Path to smalt executable;" +
                          " default: %(default)s")
    optional.add_argument("--bwa_exe", dest="bwa_exe",
                          action="store", default="bwa",
                          help="Path to BWA executable;" +
                          " default: %(default)s")
    optional.add_argument("--quast_exe", dest="quast_exe",
                          action="store", default="quast.py",
                          help="Path to quast executable; " +
                          "default: %(default)s")
    optional.add_argument("--bcftools_exe", dest="bcftools_exe",
                          action="store", default="bcftools",
                          help="Path to bcftools executable; " +
                          "default: %(default)s")
    optional.add_argument("-c", "--cores", dest='cores', action="store",
                          default=None, type=int,
                          help="cores to be used" +
                          "; default: %(default)s")
    optional.add_argument("-t", "--threads", dest='threads',
                          action="store",
                          default=1, type=int,
                          choices=[1, 2, 4],
                          help="if your cores are hyperthreaded, set number" +
                          " threads to the number of threads per processer." +
                          "If unsure, see 'cat /proc/cpuinfo' under 'cpu " +
                          "cores', or 'lscpu' under 'Thread(s) per core'." +
                          ": %(default)s")
    optional.add_argument("-m", "--memory", dest='memory', action="store",
                          default=8, type=int,
                          help="system memory available" +
                          "; default: %(default)s")
    optional.add_argument('--version', action='version',
                          version='riboSeed {version}'.format(
                              version=__version__))
    # # had to make this explicitly to call it a faux optional arg
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    args = parser.parse_args(sys.argv[2:])
    return args


def extract_node_len_cov_rc(node_name):
    rc = False
    if node_name.endswith("'"):
        rc = True
        node_name = node_name[0:-1]
    p = re.compile(r'EDGE_(?P<node>\d*?)_length_(?P<length>\d*?)_cov_(?P<cov>[\d|\.]*)')
    m = p.search(node_name)
    return (m.group("node"), int(m.group("length")), float(m.group("cov")), rc)


def make_Node(name):
    node_name, length, cov, rc = extract_node_len_cov_rc(name)
    new_node = FastgNode(
        name=node_name,
        length=length,
        cov=cov,
        reverse_complimented = rc
    )
    return new_node

def make_Neighbors(neighbor):
    node_name, length, cov, rc = extract_node_len_cov_rc(neighbor)
    new_neigh = FastgNodeNeighbor(
        name=node_name,
        reverse_complimented=rc
    )
    return new_neigh


def parse_fastg(f):
    """parse the headers in a fastg file and return a list of Node objects
    """
    node_neighs = []
    with open(f, "r") as inf:
        for line in inf:
            if line.startswith(">"):
                colons = sum([1 for x in line if x == ":" ])
                if colons > 1:
                    sys.stderr.write("multiple ':'s found in line, and can only " +
                                     "be used to separate nodes from neighbor " +
                                     "list\n")
                elif colons == 0:
                    # orphaned node or terminal node
                    # sys.stderr.write("Header does not contain a colon!\n")
                    node_neighs.append([line.strip().split(";")[0], None])
                else:
                    node_neighs.append(line.strip().split(":"))

    node_list = []
    for node, neighs in node_neighs:
        new_node = make_Node(node)
        if neighs is None:
            new_node.neighbor_list = None
        else:
            new_node.neighbor_list = [make_Neighbors(x) for x in neighs.split(",")]
        node_list.append(new_node)

    return node_list


def pathfind(node_list, top_parent, parent, prev_path, prev_length,
             path_list, thresh=1000, found_exit=False, verbose=False):
# def pathfind(node_list, parent, prev_path, prev_length, path_list, reject_direction_node, thresh=1000):
    """Returns possible exit paths from an exit node

    Given a list of all nodes, a starting node (top_parent), and information
    about any previous paths and their lengths, we recursivly trace the tree to
    find all the paths that meet a set of criteria

    1.  they originate unidirectionally from the starting node (just the
        forward or reverse compliment)
    2.  They pass through one region we deem to not be a tRNA, called exiting
    3.  Any nodes that could be within the 1000base pair threshold for flanking
        differentiation are considered
    """
    # look for both forward and rc matches
    parent_opposite_strand = [x for x in node_list if x.name == parent.name and x.reverse_complimented != parent.reverse_complimented][0]
    poss_f = [x for x in parent.neighbor_list if  x.name not in prev_path.split(":")]
    poss_rc = [x for x in parent_opposite_strand.neighbor_list if  x.name not in prev_path.split(":")]
    # this ensures we only get single direction hits from the topmost parent
    if parent == top_parent:
        possible_neighbors = poss_f
    else:
        possible_neighbors = poss_f + poss_rc
    if verbose:
        print("Prev: " + prev_path)
        print("Prev_len: " + str(prev_length))
        print("Poss:" + " ".join([x.name for x in possible_neighbors]))

    for node in possible_neighbors:
        # here we assume only one hit
        try:
            this_node = [x for x in node_list if x.name == node.name and x.reverse_complimented == node.reverse_complimented][0]
        except IndexError:
            for i in node_list:
                print("\n")
                print(i)
            print(node_list)

        if this_node.length > 250:
            found_exit = True
        if verbose:
            print("This: {0} RC: {1}".format(
                this_node.name,
                str(this_node.reverse_complimented)))
        this_length = prev_length + this_node.length
        this_path = prev_path + ":" + this_node.name
        # are we outside the zone of flanking similarity? usually 1kb?
        # and have we hit a decent stretch of sequence? one node longer than
        #   3x a tRNA, or about 270
        if this_length >= thresh and found_exit:
            if this_path in path_list:
                pass
            else:
                path_list.append(this_path)
        else:
            # further in, further in!
            # we dont capture the return list because unless this is
            # the last time, its incomplete
            pathfind(
                node_list=node_list,
                top_parent=top_parent,
                parent=this_node,
                prev_path=this_path,
                prev_length=this_length,
                path_list=path_list,
                found_exit=found_exit,
                thresh=thresh,
                verbose=verbose)
    return path_list


def make_prelim_mapping_cmds(output_root, mapping_sam, samtool_exe, spades_exe, seedGenome, k, logger):
    """ make commands to extract all reads mapping to flanking regions

    we haven't partitioned yet, but we want to do a pre-assembly of this
    partition in order to detect possile rDNA differences from the reference

    """
    region_list = [
        "{0}:{1}-{2}".format(x.sequence_id,
                             x.global_start_coord,
                             x.global_end_coord) for x in seedGenome.loci_clusters]
    logger.debug("regions to extract for the prelim mapping:\n%s",
                 "\n".join([x for x in region]))
    cmds = []
    # make a directory to contain results of this analysis
    this_dir = os.path.join(output_root, "prelim")
    os.makedirs(this_dir)
    output_sam = os.path.join(this_dir, "rDNA_reads.sam")
    output_f = os.path.join(this_dir, "rDNA_reads_1.fastq")
    output_r = os.path.join(this_dir, "rDNA_reads_2.fastq")
    output_s = os.path.join(this_dir, "rDNA_reads_s.fastq")
    # make a smatools command to extract all to a single fastq
    samtools_cmd = "{0} view -h {1} {2} > {3}".format(
        samtools_exe,
        mapping_sam,
        " ".join(regions_list),
        output_sam)
    samfastq_cmd = "{0} fastq {1} -1 {2} -2 {3} -s {4}".format(
        samtools_exe
        mapping_sam,
        output_f,
        output_r,
        output_s)
    # ensure we dont have any empty files; construct library design
    lib_count = 1
    read_libraries = ""
    for f in [output_f, output_r, output_s]:
        if os.path.getsize(f) > 0:
            read_libraries = read_libraries + " --s" + str(lib_count) + " " + f
            lib_count = lib_count + 1

    # make spades commands to run the assembly with a given k
    spades_cmd = str(
        "{0} --only-assembler --cov-cutoff off --sc --careful -k=21,33 " +
        "{1} -o {2}"
    ).format(
        spades_exe,
        read_libraries
        os.path.join(this_dir, "assembly"))
    # run the network traverse scheme

    # warn if needed


def main(args, logger=None):
    output_root = os.path.abspath(os.path.expanduser(args.output))
    try:
        os.makedirs(output_root, exist_ok=False)
    except OSError:
        print("Output directory already exists; exiting...")
        sys.exit(1)
    t0 = time.time()
    log_path = os.path.join(output_root, "riboSeed.log")
    if logger is None:
        logger = set_up_logging(verbosity=args.verbosity,
                                outfile=log_path,
                                name=__name__)

    logger.info("Usage:\n%s\n", " ".join([x for x in sys.argv]))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("%s: %s", k, str(v))

    if args.assembly_graph is None:
        if args.reference is None and args.sam is None:
            logger.error("No assembly graph provided, we must create one; " +
                         "a SAM file and a reference FASTA is required")
            sys.exit(1)
        # create a assembly graph
        assembly_graph = make_prelim_mapping_cmds(
            output_root, mapping_sam, exes, seedGenome, k, logger)

    nodes = parse_fastg(f=f)
    start = [x for x in nodes if x.name == "2" and  x.reverse_complimented ][0]
    end = [x for x in nodes if x.name == "12" and not x.reverse_complimented ][0]
    print(start)
    print(start.neighbor_list[1])
    s = pathfind(
        node_list=nodes,
        top_parent=start,
        parent=start,
        prev_path=start.name,
        prev_length=0,
        thresh=1000,
        path_list=[],
        found_exit=False)
    e = pathfind(
        node_list=nodes,
        top_parent=end,
        parent=end,
        prev_path=end.name,
        prev_length=0,
        thresh=1000,
        path_list=[],
        found_exit=False)
    print(s)
    print(e)
