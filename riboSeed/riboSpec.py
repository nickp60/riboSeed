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
DEBUG = True
PLOT = False
PLOT = True
import sys
import os
import re
import subprocess
import argparse
import multiprocessing
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot, patches
import networkx as nx

from .shared_methods import set_up_logging, make_barrnap_cmd


class FastgNode(object):
    # newid = itertools.count()

    def __init__(self, name=None, length=None, cov=None,
                 reverse_complimented=None, neighbor_list=None):
        # self.index = next(FastgNode.newid)
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
    def __init__(self, name=None, reverse_complimented=None):
        self.name = name
        self.reverse_complimented = reverse_complimented

    def __str__(self):
        return "NeighborNode: {0}\nReverse_Complimented?: {1}".format(
            self.name,
            "Yes" if self.reverse_complimented else "No"
        )



partial_list = [["EDGE_128_length_64_cov_224.111':EDGE_130_length_113_cov_395.017;", 'barrnap:0.7', 'rRNA', '2', '62', '0.00025', '+', '.', 'Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 51 percent of the 5S ribosomal RNA'], ["EDGE_128_length_64_cov_224.111:EDGE_246_length_232_cov_51.2034,EDGE_326_length_102_cov_132.553';", 'barrnap:0.7', 'rRNA', '3', '63', '0.00025', '-', '.', 'Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 51 percent of the 5S ribosomal RNA'], ["EDGE_129_length_111_cov_224.071':EDGE_130_length_113_cov_395.017;", 'barrnap:0.7', 'rRNA', '49', '109', '8.3e-05', '+', '.', 'Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 51 percent of the 5S ribosomal RNA'], ["EDGE_129_length_111_cov_224.071:EDGE_363_length_3027_cov_329.236';", 'barrnap:0.7', 'rRNA', '3', '63', '8.3e-05', '-', '.', 'Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 51 percent of the 5S ribosomal RNA'], ["EDGE_130_length_113_cov_395.017':EDGE_128_length_64_cov_224.111,EDGE_129_length_111_cov_224.071;", 'barrnap:0.7', 'rRNA', '11', '112', '1.3e-09', '-', '.', 'Name=5S_rRNA;product=5S ribosomal RNA'], ['EDGE_130_length_113_cov_395.017:EDGE_332_length_56_cov_246,EDGE_333_length_99_cov_156.886;', 'barrnap:0.7', 'rRNA', '2', '103', '1.3e-09', '+', '.', 'Name=5S_rRNA;product=5S ribosomal RNA'], ["EDGE_245_length_1702_cov_344.636':EDGE_136_length_238_cov_185.197',EDGE_244_length_141_cov_128.326;", 'barrnap:0.7', 'rRNA', '16', '1553', '0', '-', '.', 'Name=16S_rRNA;product=16S ribosomal RNA'], ['EDGE_245_length_1702_cov_344.636:EDGE_289_length_61_cov_160.167,EDGE_290_length_90_cov_217;', 'barrnap:0.7', 'rRNA', '150', '1687', '0', '+', '.', 'Name=16S_rRNA;product=16S ribosomal RNA'], ["EDGE_246_length_232_cov_51.2034':EDGE_128_length_64_cov_224.111';", 'barrnap:0.7', 'rRNA', '179', '232', '0.0033', '+', '.', 'Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 45 percent of the 5S ribosomal RNA'], ["EDGE_246_length_232_cov_51.2034:EDGE_332_length_56_cov_246';", 'barrnap:0.7', 'rRNA', '1', '54', '0.0033', '-', '.', 'Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 45 percent of the 5S ribosomal RNA'], ["EDGE_326_length_102_cov_132.553':EDGE_363_length_3027_cov_329.236';", 'barrnap:0.7', 'rRNA', '1', '54', '0.0045', '-', '.', 'Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 45 percent of the 5S ribosomal RNA'], ["EDGE_326_length_102_cov_132.553:EDGE_128_length_64_cov_224.111';", 'barrnap:0.7', 'rRNA', '49', '102', '0.0045', '+', '.', 'Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 45 percent of the 5S ribosomal RNA'], ["EDGE_363_length_3027_cov_329.236':EDGE_362_length_57_cov_268.5,EDGE_398_length_98_cov_95.6977';", 'barrnap:0.7', 'rRNA', '105', '3005', '0', '-', '.', 'Name=23S_rRNA;product=23S ribosomal RNA'], ["EDGE_363_length_3027_cov_329.236:EDGE_129_length_111_cov_224.071',EDGE_326_length_102_cov_132.553;", 'barrnap:0.7', 'rRNA', '23', '2923', '0', '+', '.', 'Name=23S_rRNA;product=23S ribosomal RNA']]

strict_list = [["EDGE_130_length_113_cov_395.017':EDGE_128_length_64_cov_224.111,EDGE_129_length_111_cov_224.071;", 'barrnap:0.7', 'rRNA', '11', '112', '1.3e-09', '-', '.', 'Name=5S_rRNA;product=5S ribosomal RNA'], ['EDGE_130_length_113_cov_395.017:EDGE_332_length_56_cov_246,EDGE_333_length_99_cov_156.886;', 'barrnap:0.7', 'rRNA', '2', '103', '1.3e-09', '+', '.', 'Name=5S_rRNA;product=5S ribosomal RNA'], ["EDGE_245_length_1702_cov_344.636':EDGE_136_length_238_cov_185.197',EDGE_244_length_141_cov_128.326;", 'barrnap:0.7', 'rRNA', '16', '1553', '0', '-', '.', 'Name=16S_rRNA;product=16S ribosomal RNA'], ['EDGE_245_length_1702_cov_344.636:EDGE_289_length_61_cov_160.167,EDGE_290_length_90_cov_217;', 'barrnap:0.7', 'rRNA', '150', '1687', '0', '+', '.', 'Name=16S_rRNA;product=16S ribosomal RNA'], ["EDGE_363_length_3027_cov_329.236':EDGE_362_length_57_cov_268.5,EDGE_398_length_98_cov_95.6977';", 'barrnap:0.7', 'rRNA', '105', '3005', '0', '-', '.', 'Name=23S_rRNA;product=23S ribosomal RNA'], ["EDGE_363_length_3027_cov_329.236:EDGE_129_length_111_cov_224.071',EDGE_326_length_102_cov_132.553;", 'barrnap:0.7', 'rRNA', '23', '2923', '0', '+', '.', 'Name=23S_rRNA;product=23S ribosomal RNA']]


def get_args():  # pragma: no cover
    """
    """
    parser = argparse.ArgumentParser(
        prog="ribo spec",
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
    optional.add_argument("-n", "--max_nodes", dest='max_nodes', action="store",
                          help="max number of nodes considered around rDNAs ",
                          type=int, default=15)
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
        name=int(node_name),
        length=length,
        cov=cov,
        reverse_complimented=rc
    )
    return new_node


def make_Neighbors(neighbor):
    node_name, length, cov, rc = extract_node_len_cov_rc(neighbor)
    new_neigh = FastgNode(
        name=int(node_name),
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
            new_node.neighbor_list = []
        else:
            new_node.neighbor_list = [make_Neighbors(x) for x in neighs.split(",")]
        node_list.append(new_node)
    return node_list


def make_adjacency_matrix(g):
    """ Make an adjacency matrix from a dict of node: [neighbors] pairs
    """
    keys=sorted(g.keys())
    size=len(keys)
    M = [ [0]*size for i in range(size) ]

    """
    for a, row in g.items() iterates over the key:value entries in dictionary, and for b in row iterates over the values. If we used (a,b), this would have given us all the pairs.

    (keys.index(a), keys.index(b)) But we need the index to assign to the corresponding matrix entry,

    keys=sorted(g.keys()) that's why we extracted and sorted the keys.

    for a,b in... getting the index entries and assigning value 1 or 2 based on diagonal element or not.

    M = [ [0]*size for ... matrix cannot be used before initialization.
    """
    for a,b in [(keys.index(a), keys.index(b)) for a, row in g.items() for b in row]:
        M[a][b] = 2 if (a==b) else 1
    return M


def draw_adjacency_matrix(G, node_order=None, partitions=[], colors=[], outdir=None):
    """
    - G is an adjacency matrix
    - node_order (optional) is a list of nodes, where each node in G
          appears exactly once
    - partitions is a list of node lists, where each node in G appears
          in exactly one node list
    - colors is a list of strings indicating what color each
          partition should be
    If partitions is specified, the same number of colors needs to be
    specified.
    """
    #Plot adjacency matrix in toned-down black and white
    fig = pyplot.figure(figsize=(4, 4)) # in inches
    pyplot.imshow(G,
                  cmap="Greys",
                  interpolation="none")

    # The rest is just if you have sorted nodes by a partition and want to
    # highlight the module boundaries
    assert len(partitions) == len(colors)
    ax = pyplot.gca()
    for partition, color in zip(partitions, colors):
        current_idx = 0
        for module in partition:
            ax.add_patch(patches.Rectangle((current_idx, current_idx),
                                          len(module), # Width
                                          len(module), # Height
                                          facecolor="none",
                                          edgecolor=color,
                                          linewidth="1"))
            current_idx += len(module)
    fig.savefig(os.path.join(outdir, "test.pdf"))


def alt_parse_fastg(f):
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
                    # jk I couldnt care less about these nodes.
                    # node_neighs.append([line.strip().split(";")[0], None])
                    node_neighs.append([line.strip()[1:-1], []])
                    pass
                else:
                    # loose the '>' at the beginning and the ';' at the end
                    node, neigh = line.strip()[1:-1].split(":")
                    node_neighs.append([node, neigh.split(",")])
    ## these should be the same length
    # print(len([x[0] for x in node_neighs]))
    # print(len(set([x[0] for x in node_neighs])))
    g = {k: v for k, v in node_neighs}
    # print([extract_node_len_cov_rc(name[0]) for name in node_neighs])
    M = make_adjacency_matrix(g)

    node_list = []
    # make objects for each node and neighbor
    for node, neighs in node_neighs:
        new_node = make_Node(node)
        if neighs is None:
            new_node.neighbor_list = []
        else:
            new_node.neighbor_list = [make_Node(x) for x in neighs]
        node_list.append(new_node)

    # make the networkx object
    DG = nx.DiGraph()
    for N in node_list:
        DG.add_node(N.name, cov=N.cov, length=N.length)
    for N in node_list:
        for neigh in N.neighbor_list:
            # if not N.reverse_complimented:
            DG.add_weighted_edges_from([(N.name, neigh.name, neigh.length)])
            # DG.add_edge(N.name, neigh.name)
            # else:
            #     DG.add_weighted_edges_from([(neigh.name, N.name, N.length)])
            # DG.add_edge(neigh.name, N.name)
    return (node_list, M, DG)


def pathfind(node_list, top_parent, parent, prev_path, prev_length,
             path_list, thresh=1000, ignored_nodes=[], found_exit=False,
             verbose=False):
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
              str(x.name) not in prev_path.split(":") and
              str(x.name) not in ignored_nodes]
    poss_rc = [x for x in parent_opposite_strand.neighbor_list if
               str(x.name) not in prev_path.split(":") and
               str(x.name) not in ignored_nodes]
    # this ensures we only get single direction hits from the topmost parent
    if parent == top_parent:
        possible_neighbors = poss_f
    else:
        possible_neighbors = poss_f + poss_rc
    if verbose:
        print("Prev: " + prev_path)
        print("Prev_len: " + str(prev_length))
        print("Poss:" + " ".join([str(x.name) for x in possible_neighbors]))

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
        this_path = prev_path + ":" + str(this_node.name)
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



def plot_G(
        G,
        nodes5,
        nodes16,
        nodes23,
        outpath, outpath2):
    fig = pyplot.figure(figsize=(10, 10)) # in inches
    pos = nx.layout.spring_layout(G, iterations=5)
    node_sizes_raw = [h['length'] for g, h in G.nodes.data()]
    maxsize = max(node_sizes_raw)
    minsize = min(node_sizes_raw)
    sizediff = maxsize - minsize
    node_sizes = [10 + (30 * ((x - minsize)/sizediff)) for x in node_sizes_raw]
    N = G.number_of_nodes()
    M = G.number_of_edges()
    _G = nx.line_graph(G)
    _N = _G.number_of_nodes()
    _M = _G.number_of_edges()
    print(N)
    print(M)
    print(_N)
    print(_M)
    edge_colors = range(2, M + 2)
    edge_alphas = [(5 + i) / (M + 4) for i in range(M)]
    # print(G.nodes.data())
    # print([h for g, h in G.nodes.data()])

    node_colors = ["lightgrey" for x in range(N)]
    for i, (g, h) in enumerate(G.nodes.data()):
        if g in nodes5:
            node_colors[i] = "red"
        if g in nodes16:
            node_colors[i] = "blue"
        if g in nodes23:
            node_colors[i] = "green"
    _edge_colors = ["lightgrey" for x in range(_M)]
    for i, (to, frm, vals) in enumerate(_G.edges.data()):
        x = [element for tupl in (to, frm) for element in tupl]
        if len(set(x).intersection(nodes5)) > 0:
            _edge_colors[i] = "red"
        if len(set(x).intersection(nodes16)) > 0:
            _edge_colors[i] = "blue"
        if len(set(x).intersection(nodes23)) > 0:
            _edge_colors[i] = "green"
        # node_colors
    nx.draw(G, with_labels=True, linewidths=0, node_size=node_sizes,
            alpha=0.7,  font_size=2, arrows=True,
            node_color=node_colors, edge_color="darkgrey", width=.2)
    # edges = nx.draw_networkx_nodes(G, pos,
    #                                linewidths=0,
    #                                with_labels=True,
    #                                node_size=node_sizes,
    #                                node_color=node_colors)
    # nodes = nx.draw_networkx_edges(G, pos, node_size=node_sizes, arrowstyle='->',
    #                                arrowsize=5, edge_color="darkgrey",
    #                                with_labels=True,
    #                                edge_cmap=plt.cm.Blues, width=.3)
    # set alpha value for each edge
    # for i in range(M):
    #     edges[i].set_alpha(edge_alphas[i])
    ax = plt.gca()
    ax.set_axis_off()
    fig.savefig(outpath)


    fig = pyplot.figure(figsize=(10, 10)) # in inches
    nx.draw(nx.line_graph(G), with_labels=True, linewidths=0, node_size=node_sizes,
            alpha=0.7,  font_size=2, arrows=False,
            node_color="black", edge_color=_edge_colors, width=.2)
    fig.savefig(outpath2)

def return_if_interior_node(node, lengths, cutoff, verbose=False):
    """ return tuple of (valid_node, IS_BORDER)

    See docstring for neighborhood_by_length for greater detail

    """
    if lengths[node] > cutoff:
        return (None, None)
    else:
        if verbose:
            print("found interior node: %s" %node)
        return (node, False)


def return_if_border_node(node, paths, interior_nodes):
    """ return tuple of (valid_node, IS_BORDER)

    See docstring for neighborhood_by_length for greater detail

    """
    IS_BORDER = False
    for p in paths[node]:
        if p in interior_nodes:
            IS_BORDER = True
    if IS_BORDER:
        print("found border node %s " % node)
        return (node, True)
    else:
        return (None, True)


def neighborhood_by_n(G, node, n):
    path_lengths = nx.single_source_dijkstra_path_length(G, node)
    return [node for node, length in path_lengths.items()
            if length == n]


def neighborhood_by_length(G, source, cutoff=20000, ignored_nodes=[]):
    """
    I needed a way to see if a given node was within a certain distance from a source node by the shortest path.  This could be done with the dijkstra_predecessor_and_distance function, but that rejects any path >= the cutoff, whereas I need to retain the ones where the cutoff occurs within the path too.  So this takes a list of ALL the paths and their lengths (from  the dictionaries returned by networkx's dijkstra_predecessor_and_distance used without a cutoff), and then recursivly iterates through the network.

    for each node in the paths returned by dijkstra_predecessor_and_distance,
    we run through the find_inclusive_paths_within_cutoff method, calls nodes
    as either being interior (ie, within the cutoff), or on the border
    (ie, the edge included the cutoff)

    """
    # print(nx.single_source_dijkstra_path_length(G, source))
    # print(path_nodes)
    # print(path_lengths)
    interior_nodes = [source]
    border_nodes = []
    nodesDict = dict(G.nodes(data=True))
    # sys.exit()
    paths_dict = nx.single_source_dijkstra(G, source)[1]
    # print(paths_dict)
    for target, path_to in paths_dict.items():
        path_len = 0
        if len(set(path_to).intersection(set(ignored_nodes))) > 0:
            continue
        for i, node in enumerate(path_to):
            if i > 0:
                path_len = path_len + nodesDict[node]['length']
                if path_len > cutoff:
                    border_nodes.append(node)
                    break
                elif path_len < cutoff:
                    interior_nodes.append(node)
                else:
                    pass
    #     included_node, is_border = return_if_interior_node(
    #         node=n,
    #         lengths=path_lengths,
    #         cutoff=cutoff)
    #     # included_node will be None if we have already dealt with it
    #     if included_node is not None:
    #             interior_nodes.append(included_node)

    # for n, p in path_nodes.items():
    #     if n == source or n in interior_nodes:
    #         continue
    #     included_node, is_border = return_if_border_node(
    #         node=n,
    #         paths=path_nodes,
    #         interior_nodes=interior_nodes)
    #     # included_node will be None if we have already dealt with it
    #     if included_node is not None:
    #             border_nodes.append(included_node)
    # print(set(interior_nodes))
    # print(set(border_nodes))
    return (set(interior_nodes), set(border_nodes))


def make_gff_list(gffpath):
    gff_list = []
    with open(gffpath, 'r') as g:
        for idx, line in enumerate(g):
            if idx == 0:
                #gets rid of header line
                pass
            else:
                gff_list.append(line.strip().split("\t"))
    return gff_list


def find_rRNA_from_gffs(gff_list, partial=False):
    """
    this is a bit convoluted
    for gff lines where the product has the nane of the rDNA (16S, 23S, etc),
      we extract the node name (usually a number), which is the first thing
      returned from the extract_node_len_cov function
    """
    nodes16, nodes23, nodes5, = [], [], []
    print(gff_list)
    for x in gff_list:
        if ("partial" in x[8] and not partial):
            continue
        if "16S" in x[8]:
            nodes16.append(int(extract_node_len_cov_rc(x[0])[0]))
        if "23S" in x[8]:
            nodes23.append(int(extract_node_len_cov_rc(x[0])[0]))
        if "5S" in x[8]:
            nodes5.append(int(extract_node_len_cov_rc(x[0])[0]))
    nodes16, nodes23, nodes5 = list(set(nodes16)), list(set(nodes23)), list(set(nodes5))
    print("16s nodes:")
    print(nodes16)
    print("23s nodes:")
    print(nodes23)
    print("5s nodes:")
    print(nodes5)
    return(nodes16, nodes23, nodes5)



# def coalesce(G,node1,node2):
#     """Performs Briggs coalescing. Takes in the graph and two nodes.
#     Returns 1 if unable to coalesce, 0 otherwise.
#     https://stackoverflow.com/questions/17483022
#     """
#     if node1 in G.neighbors(node2) or node2 in G.neighbors(node1):
#         print "Cannot coalesce. Node",node1,"and node",node2,"share an edge"
#         return 1
#     elif G.degree(node1)+G.degree(node2) >= k:
#         print "Cannot coalesce. Combined degree of",node1,"and",node2,"\
# is",G.degree(node1)+G.degree(node2),"which is too high for k =",k
#         return 1
#     else:
#         newedge = []
#         for i in range(len(G.neighbors(node2))):
#             newedge.append((node1 , G.neighbors(node2)[i]))
#         G.add_edges_from(newedge)
#         G.remove_node(node2)
#         nx.relabel_nodes(G, {node1:node1+node2},copy=False)
#     return 0


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
    nodes, M, G = alt_parse_fastg(f=args.assembly_graph)

    draw_adjacency_matrix(M, node_order=None, partitions=[], colors=[], outdir=args.output)
    #  write out the adjacency matrix
    with open(os.path.join(args.output, "tab.txt"), "w") as o:
        for line in M:
            o.write("\t".join([str(x) for x in line]) + "\n")
    # alt nodes
    dic = {int(x.name): [int(y.name) for y in x.neighbor_list] for x in nodes}


    # run barrnap to find our rDNAs; first time trhough is to find "full" -
    # length genes, and the second run with the relaxed thresholds is for
    # finding partial genes
    barrnap_gff = os.path.join(output_root, "strict_barrnapped.gff")
    barrnap_gff_partial = os.path.join(output_root, "partial_barrnapped.gff")
    barrnap_cmd = make_barrnap_cmd(
        infasta=args.assembly_graph,
        outgff=barrnap_gff,
        exe=args.barrnap_exe,
        threads=args.cores,
        thresh=0.8,
        evalue=1e-06,
        kingdom="bac")
    barrnap_cmd_partial = make_barrnap_cmd(
        infasta=args.assembly_graph,
        outgff=barrnap_gff_partial,
        exe=args.barrnap_exe,
        threads=args.cores,
        thresh=0.1,
        evalue=1,
        kingdom="bac")
    if DEBUG:
        gff_list = strict_list
        gff_list_partial = partial_list
    else:
        for cmd in [barrnap_cmd, barrnap_cmd_partial]:
            logger.info("running barrnap cmd: %s", cmd)
            subprocess.run(cmd,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
        # determine which ones are our 16s, 23s, and 5s nodes
        gff_list = make_gff_list(barrnap_gff)
        logger.debug(gff_list)
        gff_list_partial = make_gff_list(barrnap_gff_partial)
        logger.debug(gff_list_partial)

    solid16, solid23, solid5 = find_rRNA_from_gffs(gff_list, partial=False)
    partial16, partial23, partial5 = find_rRNA_from_gffs(gff_list_partial, partial=True)
    partial16 = [x for x in partial16 if x not in solid16]
    partial23 = [x for x in partial23 if x not in solid23]
    partial5 = [x for x in partial5 if x not in solid5]
    rrnas = {
        "16S": {
            "partial": partial16,
            "solid": solid16
        },
        "23S": {
            "partial": partial23,
            "solid": solid23
        },
        "5S": {
            "partial": partial5,
            "solid": solid5
        }
    }
    print(rrnas)
    if len(solid16) > 1:
        print("more than one full 16S contig found")
    if len(solid23) > 1:
        print("more than one full 23S contig found")





    oldG = deepcopy(G)


    ################
    # collapse nodes where partial loci neighbor full-length
    collapsed = []
    print(rrnas)
    for k, vals in rrnas.items():
        print("checking for collapsable %s nodes" %k)
        these_collapsed = []
        if len(vals["partial"]) == 0:
            continue
        # check if partial node neighbors true node
        # print(vals)
        for part in vals["partial"]:
            # print("partial: %i" %part)
            for solid in vals["solid"]:
                if part in G.neighbors(solid):
                    # print(G.edges(solid))
                    for d in G.edges(part):
                        if d[1] != solid and d[1] not in G.neighbors(solid):
                            print(d)
                            partial_to_next_weight = G.get_edge_data(d[0],d[1])["weight"]
                            print(partial_to_next_weight)
                            solid_to_partial_weight = G.get_edge_data(solid, d[0])["weight"]
                            print(solid_to_partial_weight)
                            G.add_edge(solid , d[1], weight = solid_to_partial_weight + partial_to_next_weight)
                    G.remove_node(part)
                    these_collapsed.append(part)
        # print(G.get_edge_data(solid, d[1]))
        print("removed %i nodes:" % len(these_collapsed))
        print(these_collapsed)
        collapsed.extend(these_collapsed)
    if PLOT:
        plot_G(
            G,
            solid5,
            solid16,
            solid23,
            outpath=os.path.join(args.output, "test_post_collapse_G.pdf"),
            outpath2=os.path.join(args.output, "test_post_collapse_G_linegraph.pdf"),
        )

    ########  Reduce this graph to all nodes within 20kb if a 16S region
    interior_nodes = []
    border_nodes = []
    for i in solid16:
        # print(i)
        # print([x for x in G.neighbors(i)])
        # first = max(x for x in G.neighbors(i))
        # len_16S = G.get_edge_data(i, first)
        # print([G.get_edge_data(i, x) for x in G.neighbors(i)])
        # print(len_16S)
        interior, border = neighborhood_by_length(G, i, cutoff=1000, ignored_nodes=solid23)
        interior_nodes.extend(interior)
        border_nodes.extend(border)
    for i in solid23:
        interior, border = neighborhood_by_length(G, i, cutoff=1000, ignored_nodes=solid16)
        interior_nodes.extend(interior)
        border_nodes.extend(border)
    valid_nodes = [x for y in [interior_nodes, border_nodes] for x in y]

    print(len(set(valid_nodes)))
    bad_nodes = set(G.nodes).symmetric_difference(set(valid_nodes))
    print(len(G.nodes))
    print(len(bad_nodes))
    for node in bad_nodes:
        G.remove_node(node)
    if PLOT:
        plot_G(
            G,
            solid5,
            solid16,
            solid23,
            outpath=os.path.join(args.output, "test_post_reduction_G.pdf"),
            outpath2=os.path.join(args.output, "test_post_reduction_G_linegraph.pdf"),
        )


    #######   Collapse shtuff between the 16S and 23S, if its less than say 2kb
    # for k, vals in rrnas.items():
    #     print("checking for collapsable %s nodes" %k)
    #     these_collapsed = []
    #     if len(vals["partial"]) == 0:
    #         continue
    #     # chekc if partial node neighbors true node
    #     # print(vals)
    #     for part in vals["partial"]:
    #         # print("partial: %i" %part)
    connector_paths = []
    for node16 in rrnas["16S"]["solid"]:
        for node23 in rrnas["23S"]["solid"]:
            connector_paths.extend(
                nx.all_simple_paths(G, node16, node23, cutoff=500)
            )
    print("number of connector paths: %i" %len(connector_paths))
    # count the paths going out from the 16S
    out_paths = []
    print([G.out_degree(g) for g in G.nodes])
    tips = [node for node in G.nodes() if G.out_degree(node) == 1]
    print("tips:")
    print(tips)
    for node16 in rrnas["16S"]["solid"]:
        for tip in tips:
            out_paths.extend(nx.all_simple_paths(G, node16, tip))
    print("number of out paths: %i" %len(out_paths))
    for path in connector_paths:
        l = len(path)
        for opath in out_paths:
            if path == opath[0:l]:
                out_paths.remove(opath)
    print("number of filtered out paths: %i" %len(out_paths))
    print("23S outpaths")
    print("number of connector paths: %i" %len(connector_paths))
    # count the paths going out from the 16S
    out_paths_23 = []
    for node23 in rrnas["23S"]["solid"]:
        for tip in tips:
            out_paths_23.extend(nx.all_simple_paths(G, node23, tip))
    print("number of 23S out paths: %i" %len(out_paths_23))
    for path in connector_paths:
        l = len(path)
        rpath = list(reversed(path))
        for opath in out_paths_23:
            if 245 in opath:
            # if rpath == opath[-l: ]:
                out_paths_23.remove(opath)
    print("number of 23S filtered out paths: %i" %len(out_paths_23))
    for i in out_paths_23:
        print(i)


    ###
    # plot_G(
    #     oldG,
    #     solid5,
    #     solid16,
    #     solid23,
    #     outpath=os.path.join(args.output, "test_oldG.pdf"),
    #     outpath2=os.path.join(args.output, "test_oldG_linegraph.pdf"),
    # )
    ########
    sys.exit()
    if len(rrnas["16S"]["solid"]) >  1:
        logger.error("it appears that there are distinct 16S rDNAs in the " +
                     "assenbly graph; this tracing algorithm is not the best" +
                     "option.  Please review the graph manually to determine" +
                     "probable number of rDNAs")
        raise ValueError
    elif len(rrnas["16S"]["solid"]) < 1:
        logger.error("Barrnap failed to detect any full 16S in the assembly graph")
        raise ValueError

    if len(rrnas["23S"]["solid"]) > 1:
        logger.error("it appears that there are distinct 23S rDNAs in the " +
                     "assenbly graph; this tracing algorithm is not the best" +
                     "option.  Please review the graph manually to determine" +
                     "probable number of rDNAs")
        raise ValueError
    elif len(rrnas["23S"]["solid"]) < 1:
        logger.error("Barrnap failed to detect any 23S in the assembly graph")
        raise ValueError

    # now, we set out start nodes to be 16S and end to be 23S
    start = [x for x in nodes if x.name == node16 and  x.reverse_complimented ][0]
    logger.info("16S node: %s", start)
    end = [x for x in nodes if x.name == node23 and not x.reverse_complimented ][0]
    logger.info("23S node: %s", end)
    s = pathfind(
        node_list=nodes,
        top_parent=start,
        parent=start,
        prev_path=str(start.name),
        prev_length=0,
        thresh=1000,
        ignored_nodes=[node23],
        path_list=[],
        found_exit=False)
    e = pathfind(
        node_list=nodes,
        top_parent=end,
        parent=end,
        prev_path=str(end.name),
        prev_length=0,
        thresh=2000,
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
