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


def draw_adjacency_matrix(G, node_order=None, partitions=[], colors=[], outdir=None):
    """
    - G is a dictionary that can be comnverted to a network graph
    - node_order (optional) is a list of nodes, where each node in G
          appears exactly once
    - partitions is a list of node lists, where each node in G appears
          in exactly one node list
    - colors is a list of strings indicating what color each
          partition should be
    If partitions is specified, the same number of colors needs to be
    specified.
    """
    # G = nx.Graph(g)
    adjacency_matrix = G
    # adjacency_matrix = nx.to_numpy_matrix(G, dtype=np.bool, nodelist=node_order)

    #Plot adjacency matrix in toned-down black and white
    fig = pyplot.figure(figsize=(4, 4)) # in inches
    pyplot.imshow(adjacency_matrix,
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


# TODO replace the matrix with a adjacency matrix
# https://stackoverflow.com/questions/37353759/how-do-i-generate-an-adjacency-matrix-of-a-graph-from-a-dictionary-in-python
# replace pathfinding algo with a matrix traverse rather than this spitshow

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



def plot_G(G,         nodes5,
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
            node_colors[i] = "blue"
        if g in nodes16:
            node_colors[i] = "red"
        if g in nodes23:
            node_colors[i] = "green"
    _edge_colors = ["lightgrey" for x in range(_M)]
    for i, (to, frm, vals) in enumerate(_G.edges.data()):
        x = [element for tupl in (to, frm) for element in tupl]
        if len(set(x).intersection(nodes5)) > 0:
            _edge_colors[i] = "blue"
        if len(set(x).intersection(nodes16)) > 0:
            _edge_colors[i] = "red"
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

def return_if_interior_node(source, node, lengths, paths, cutoff,
                            interior_nodes, border_nodes):
    """ return tuple of (valid_node, IS_BORDER)

    See docstring for neighborhood_by_length for greater detail

    """
    if lengths[node] > cutoff:
        return (None, None)
    else:
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


def neighborhood_by_length(G, source, cutoff=20000):
    """
    I needed a way to see if a given node was within a certain distance from a source node by the shortest path.  This could be done with the dijkstra_predecessor_and_distance function, but that rejects any path >= the cutoff, whereas I need to retain the ones where the cutoff occurs within the path too.  So this takes a list of ALL the paths and their lengths (from  the dictionaries returned by networkx's dijkstra_predecessor_and_distance used without a cutoff), and then recursivly iterates through the network.

    for each node in the paths returned by dijkstra_predecessor_and_distance,
    we run through the find_inclusive_paths_within_cutoff method, calls nodes
    as either being interior (ie, within the cutoff), or on the border
    (ie, the edge included the cutoff)

    """
    path_nodes, path_lengths = nx.dijkstra_predecessor_and_distance(G, source)
    # print(nx.single_source_dijkstra_path_length(G, source))
    # print(path_nodes)
    # print(path_lengths)
    interior_nodes = [source]
    border_nodes = []
    for n, p in path_nodes.items():
        if n == source:
            continue
        included_node, is_border = return_if_interior_node(
            source=source,
            node=n,
            lengths=path_lengths,
            paths=path_nodes,
            cutoff=cutoff,
            interior_nodes=interior_nodes,
            border_nodes=border_nodes)
        # included_node will be None if we have already dealt with it
        if included_node is not None:
                interior_nodes.append(included_node)

    for n, p in path_nodes.items():
        if n == source or n in interior_nodes:
            continue
        included_node, is_border = return_if_border_node(
            node=n,
            paths=path_nodes,
            interior_nodes=interior_nodes)
        # included_node will be None if we have already dealt with it
        if included_node is not None:
                border_nodes.append(included_node)

    print(interior_nodes)
    print(border_nodes)
    return (interior_nodes, border_nodes)


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
    with open(os.path.join(args.output, "tab.txt"), "w") as o:
        for line in M:
            o.write("\t".join([str(x) for x in line]) + "\n")
    # alt nodes
    dic = {int(x.name): [int(y.name) for y in x.neighbor_list] for x in nodes}
    # print(dic)


    # run barrnap to find our rDNAs
    barrnap_gff = os.path.join(output_root, "barrnapped.gff")
    barrnap_cmd = make_barrnap_cmd(
        infasta=args.assembly_graph,
        outgff=barrnap_gff,
        exe=args.barrnap_exe,
        threads=args.cores,
        thresh=0.1,
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
    nodes16, nodes23, nodes5, = [], [], []
    print(gff_list)
    for x in gff_list:
        print(x)
        if "16S" in x[8]:
            nodes16.append(int(extract_node_len_cov_rc(x[0])[0]))
        if "23S" in x[8]:
            nodes23.append(int(extract_node_len_cov_rc(x[0])[0]))
        if "5S" in x[8]:
            nodes5.append(int(extract_node_len_cov_rc(x[0])[0]))
    print(nodes16)
    print(nodes23)
    print(nodes5)

    # print(nodes5)
    # print(color_mask)
    ########  Reduce this graph
    oldG = deepcopy(G)
    interior_nodes = []
    border_nodes = []
    # max_depth = 15

    interior, border = neighborhood_by_length(G, nodes16[0], cutoff=1000)
    interior_nodes.extend(interior)
    border_nodes.extend(border)
    valid_nodes = [x for y in [interior_nodes, border_nodes] for x in y]

    print(len(set(valid_nodes)))
    bad_nodes = set(G.nodes).symmetric_difference(set(valid_nodes))
    print(len(G.nodes))
    print(len(bad_nodes))
    for node in bad_nodes:
        G.remove_node(node)
    #########
    plot_G(
        G,
        nodes5,
        nodes16,
        nodes23,
        outpath=os.path.join(args.output, "test_G.pdf"),
        outpath2=os.path.join(args.output, "test_G_linegraph.pdf"),
    )
    plot_G(
        oldG,
        nodes5,
        nodes16,
        nodes23,
        outpath=os.path.join(args.output, "test_oldG.pdf"),
        outpath2=os.path.join(args.output, "test_oldG_linegraph.pdf"),
    )
    ########
    sys.exit()
    if not nodes16.count(nodes16[0]) == len(nodes16):
        logger.error("it appears that there are distinct 16S rDNAs in the " +
                     "assenbly graph; this tracing algorithm is not the best" +
                     "option.  Please review the graph manually to determine" +
                     "probable number of rDNAs")
        raise ValueError
    elif len(nodes16) < 1:
        logger.error("Barrnap failed to detect any 16S in the assembly graph")
        raise ValueError

    if len(nodes23) == 0 or not nodes23.count(nodes23[0]) == len(nodes23):
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
