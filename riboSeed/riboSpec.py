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
import subprocess
import argparse
import multiprocessing

from .shared_methods import set_up_logging, make_barrnap_cmd, test_barrnap_ok


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
        description="Given either an assembly graph or a mapping file " +
        "and reference, determine whether the number of rDNAs appears " +
        "to match the reference",
        add_help=False)  # to allow for custom help
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-g", "--assemgly_graph",
                          dest='assembly_graph',
                          action="store", default='', type=str,
                          help="fastg assembly graph from SPAdes",
                          required=False)
    optional.add_argument("-b", "--bam", dest='mapping_bam', action="store",
                          help="indexed mapping file of reads to the " +
                          "reference; required if no assembly graph is used",
                          type=str, default=None)
    optional.add_argument("-r", "--reference", dest='reference', action="store",
                          help="the same reference fasta used to generate " +
                          "mapping file; "
                          "required if no assembly graph is used",
                          type=str, default=None)
    optional.add_argument("-v", "--verbosity", dest='verbosity',
                          action="store",
                          default=2, type=int, choices=[1, 2, 3, 4, 5],
                          help="Logger writes debug to file in output dir; " +
                          "this sets verbosity level sent to stderr. " +
                          " 1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
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
    optional.add_argument("--barrnap_exe", dest="barrnap_exe",
                          action="store", default="barrnap",
                          help="Path to barrnap executable;" +
                          " default: %(default)s")
    optional.add_argument("-c", "--cores", dest='cores', action="store",
                          default=None, type=int,
                          help="cores to be used" +
                          "; default: %(default)s")
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
             path_list, thresh=1000, ignored_nodes=[], found_exit=False, verbose=False):
    """Returns possible exit paths from an exit node

    Given a list of all nodes, a starting node (top_parent), and information
    about any previous paths and their lengths, we recursivly trace the tree to
    find all the paths that meet a set of criteria

    1.  they originate unidirectionally from the starting node (just the
        forward or reverse compliment)
    2.  They pass through one region we deem to not be a tRNA, called exiting
    3.  Any nodes that could be within the 1000base pair threshold for flanking
        differentiation are considered
    4.  They dont pass through ignored_nodes, which alows us to set the directionality
        (ie, paths from the 16S cant pass through a 23S)
    """
    # look for both forward and rc matches
    parent_opposite_strand = [x for x in node_list if x.name == parent.name and x.reverse_complimented != parent.reverse_complimented][0]
    poss_f = [x for x in parent.neighbor_list if
              x.name not in prev_path.split(":") and
              x.name not in ignored_nodes]
    poss_rc = [x for x in parent_opposite_strand.neighbor_list if
               x.name not in prev_path.split(":") and
               x.name not in ignored_nodes]
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


def run_prelim_mapping_cmds(output_root, mapping_sam, samtool_exe, spades_exe, seedGenome, k, logger):
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
        samtools_exe,
        mapping_sam,
        output_f,
        output_r,
        output_s)
    logger.debug("running commands to get reads mapping to any rDNA region")
    for cmd in [samtools_cmd, samfastq_cmd]:
        logger.debug(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
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
        read_libraries,
        os.path.join(this_dir, "assembly"))
    logger.debug("running commands rDNA-mapping reads")
    logger.debug(spades_cmd)
    subprocess.run(spades_cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)

def make_simple_header():
    """uyse sed to just get node name from fastg
    sed 's/^[^ ]\(.*\)[:]\(.*\).*$/>\1/' assembly_graph.fastg > renamed.fastg
    """
    pass


def main(args, logger=None):
    output_root = os.path.abspath(os.path.expanduser(args.output))
    try:
        os.makedirs(output_root, exist_ok=False)
    except OSError:
        print("Output directory already exists; exiting...")
        sys.exit(1)
    log_path = os.path.join(output_root, "riboSpec.log")
    if logger is None:
        logger = set_up_logging(verbosity=args.verbosity,
                                outfile=log_path,
                                name=__name__)

    logger.info("Usage:\n%s\n", " ".join([x for x in sys.argv]))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("%s: %s", k, str(v))
    if args.cores is None:
        args.cores = multiprocessing.cpu_count()
    if args.assembly_graph is None:
        if args.reference is None and args.sam is None:
            logger.error("No assembly graph provided, we must create one; " +
                         "a SAM file and a reference FASTA is required")
            sys.exit(1)
        # create a assembly graph
        args.assembly_graph = make_prelim_mapping_cmds()
    # make a list of node objects
    nodes = parse_fastg(f=args.assembly_graph)

    # run barrnap to find our rDNAs
    barrnap_gff = os.path.join(output_root, "barrnapped.gff")
    barrnap_cmd = make_barrnap_cmd(
        infasta=args.assembly_graph,
        outgff=barrnap_gff,
        exe=args.barrnap_exe,
        threads=args.cores,
        thresh=0.8,
        kingdom="bac")
    logger.info("running barrnap cmd: %s", barrnap_cmd)
    subprocess.run(barrnap_cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    # determine which ones are our 16s, 23s, and 5s nodes
    gff_list = []
    with open(barrnap_gff, 'r') as g:
        for idx, line in enumerate(g):
            if idx == 0:
                #gets rid of header line
                pass
            else:
                gff_list.append(line.strip().split("\t"))
    logger.debug(gff_list)
    # this is a bit convoluted
    # for gff lines where the product has the nane of the rDNA (16S, 23S, etc),
    #   we extract the node name (usually a number), which is the first thing
    #   returned from the extract_node_len_cov function
    nodes16 = [ extract_node_len_cov_rc(x[0])[0] for x in gff_list if "16S" in x[8]]
    nodes23 = [ extract_node_len_cov_rc(x[0])[0] for x in gff_list if "23S" in x[8]]
    nodes5  = [ extract_node_len_cov_rc(x[0])[0] for x in gff_list if "5S"  in x[8]]
    #
    if not nodes16.count(nodes16[0]) == len(nodes16):
        logger.error("it appears that there are distinct 16S rDNAs in the " +
                     "assenbly graph; this tracing algorithm is not the best" +
                     "option.  Please review the graph manually to determine" +
                     "probable number of rDNAs")
        raise ValueError
    elif len(nodes16) < 1:
        logger.error("Barrnap failed to detect any 16S in the assembly graph")
        raise ValueError

    if not nodes23.count(nodes23[0]) == len(nodes23):
        logger.error("it appears that there are distinct 23S rDNAs in the " +
                     "assenbly graph; this tracing algorithm is not the best" +
                     "option.  Please review the graph manually to determine" +
                     "probable number of rDNAs")
        raise ValueError
    elif len(nodes23) < 1:
        logger.error("Barrnap failed to detect any 23S in the assembly graph")
        raise ValueError
    node16 = nodes16[0]
    node23 = nodes23[0]

    # now, we set out start nodes to be 16S and end to be 23S
    start = [x for x in nodes if x.name == node16 and  x.reverse_complimented ][0]
    logger.info("16S node: %s", start)
    end = [x for x in nodes if x.name == node23 and not x.reverse_complimented ][0]
    logger.info("23S node: %s", end)
    s = pathfind(
        node_list=nodes,
        top_parent=start,
        parent=start,
        prev_path=start.name,
        prev_length=0,
        thresh=1000,
        ignored_nodes=[node23],
        path_list=[],
        found_exit=False)
    e = pathfind(
        node_list=nodes,
        top_parent=end,
        parent=end,
        prev_path=end.name,
        prev_length=0,
        thresh=1000,
        ignored_nodes=[node16],
        path_list=[],
        found_exit=False)
    logger.info("Unique possible paths exiting the 16S rDNA: \n\t{0}".format(
        "\n\t".join(s)))
    logger.info("Unique possible paths exiting the 16S rDNA: \n\t{0}".format(
        "\n\t".join(e)))

    # interpret results
    n_upstream = len(s)
    n_downstream = len(e)
    logger.info("Paths leading to rDNA operon: %i", n_upstream)
    logger.info("Paths exiting  rDNA operon: %i", n_downstream)
    if n_upstream == n_downstream:
        logger.info("This indicates that there are at least %i rDNAs present",
                    n_upstream)
    else:
        logger.info("inconclusive results")
