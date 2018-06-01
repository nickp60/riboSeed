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

"""

DEBUG = False
# DEBUG = True
# PLOT = False
# PLOT = True
import sys
import os
import re
import random
import math
import subprocess
import argparse
import multiprocessing
from copy import deepcopy


try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import pyplot, patches
    PLOT = True
except Exception as e:  # likely an ImportError, but not taking chances
    print(e)
    print("\nlooks like you have some issue with matplotlib.  " +
          "Classic matplotlib, amirite? Plotting is disabled\n")
    PLOT = False

import numpy as np
import networkx as nx


from .shared_methods import set_up_logging, make_barrnap_cmd


class FastgNode(object):
    def __init__(self, name=None, length=None, cov=None,
                 reverse_complimented=None, neighbor_list=None, raw=None):
        # self.index = next(FastgNode.newid)
        self.name = name
        self.length = length
        self.cov = cov
        self.neighbor_list = neighbor_list
        self.reverse_complimented = reverse_complimented
        self.raw = raw

    def __str__(self):
        return str("Node: {0}\nNeighbors: {1}\nLength: {2}\nCoverage: " +
                   "{3}\nReverse_Complimented?: {4}\nRaw: {5} ").format(
                       self.name,
                       "None" if self.neighbor_list is None else
                           ",".join([str(x.name) for x in self.neighbor_list]),
                       "None" if self.length is None else self.length,
                       "None" if self.cov is None else self.cov,
                       str(self.reverse_complimented),
                       "None" if self.raw is None else self.raw,
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
    requiredNamed.add_argument(
        "-o", "--output", dest='output', action="store",
        help="output directory; " +
        "default: %(default)s", default=os.getcwd(),
        type=str, required=True)
    requiredNamed.add_argument(
        "-g", "--assemgly_graph",
        dest='assembly_graph',
        action="store", default='', type=str,
        help="fastg assembly graph from SPAdes",
        required=True)
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--plot_graphs", dest='plot_graphs',
                          help="draw the network graphs ",
                          action="store_true")
    optional.add_argument("-v", "--verbosity", dest='verbosity',
                          action="store",
                          default=2, type=int, choices=[1, 2, 3, 4, 5],
                          help="Logger writes debug to file in output dir; " +
                          "this sets verbosity level sent to stderr. " +
                          " 1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
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
        reverse_complimented=rc,
        raw=name
    )
    return new_node


def make_Neighbors(neighbor):
    node_name, length, cov, rc = extract_node_len_cov_rc(neighbor)
    new_neigh = FastgNode(
        name=int(node_name),
        reverse_complimented=rc
    )
    return new_neigh


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


def plot_adjacency_matrix(G, node_order=None, partitions=[], colors=[], outpath=None):
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
    fig.savefig(outpath)


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
    return (node_list, M, DG)


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
    ax = pyplot.gca()
    ax.set_axis_off()
    fig.savefig(outpath)


    fig = pyplot.figure(figsize=(10, 10)) # in inches
    nx.draw(nx.line_graph(G), with_labels=True, linewidths=0, node_size=node_sizes,
            alpha=0.7,  font_size=2, arrows=False,
            node_color="black", edge_color=_edge_colors, width=.2)
    fig.savefig(outpath2)


def neighborhood_by_n(G, node, n):
    path_lengths = nx.single_source_dijkstra_path_length(G, node)
    return [node for node, length in path_lengths.items()
            if length == n]


def get_nodes_weight_instead_of_edges(u, v, d):
    print(u)
    print(v)
    print(d)
    # ie from  https://stackoverflow.com/questions/49136427/
    # node_u_wt = G.nodes[u].get('node_weight', 1)
    node_v_wt = G.nodes[v].get('length', 1)
    # edge_wt = d.get('weight', 1)
    return node_v_wt




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
    # changed this from dijkstra path because that was getting the shortes path by number of intermnediate nodes.  We can only weight nodes, not edges, so this was problematic
    paths_dict = {}
    # first, use single_source_dij to get all connected nodes.
    get_node_weight = lambda u,v,d: G.node[v].get('length', 0)

    targets_path = nx.single_source_dijkstra(G, source, weight=get_node_weight)[0]
    with open("tmp_paths.txt", "w") as out:
        for target, path_to in paths_dict.items():
            out.write(str(target) + "     " + " ".join([str(x) for x in path_to]) + "\n")
    sys.exit()
    if 0:
        path_len = 0
        # ignore paths that travel through the "bad places"
        if len(set(path_to).intersection(set(ignored_nodes))) > 0:
            # print("ignoring path")
            # print(path_to)                  
            continue
        for i, node in enumerate(path_to):
            if i > 0: # ignore source
                path_len = path_len + nodesDict[node]['length']
                if path_len > cutoff:
                    border_nodes.append(node)
                    break
                else:
                    interior_nodes.append(node)
        print(path_to)
        print(path_len)
    sys.exit()
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


def get_depth_of_big_nodes(G, threshold=5000):
    nodes_dict = dict(G.nodes(data=True))
    lengths = []
    depths = []
    prods = []
    stupid_normalized_list = []
    for node, data in nodes_dict.items():
        if data['length'] > threshold:
            depths.append(data['cov'])
            lengths.append(data['length'])
    for i, l in enumerate(lengths):
        prods.append(depths[i] * l)
        # this is to get stats from incorporating some notion of weight
        stupid_normalized_list.extend([depths[i] for x in range(math.ceil(l/1000))])
    # print("normalized:")
    # print(stupid_normalized_list)
    totlen = sum(lengths)
    ave = sum([x for x in prods]) /totlen
    print("Total length of nodes passing threshold: %i" % totlen)
    print("Average depth of contigs greater %i is %f" %(threshold, ave))
    quarts = []
    for quart in [.25, .5, .75]:
        q = percentile(stupid_normalized_list,quart)
        print("%i th quartile: %f" %  \
              (quart * 100, q))
        quarts.append(q)
    return(depths, ave, quarts)


def make_silly_boxplot(vals, outpath, names=None, title="", yline=None):
    """
    yline can either be an number or a tuple: (25%, 50%, 75% quartiles)
    """
    fig = pyplot.figure(figsize=(6, 6)) # in inches
    ax = pyplot.axes()
    if isinstance(vals[0], list):
        pyplot.boxplot(vals, 0, 'rs', 1)
        if names is not None:
            assert len(names) == len(vals), \
                "number of names does not equal number of sets!"
            ax.set_xticklabels(names)

        for i, data in enumerate(vals) :
            for d in data:
                pyplot.scatter(x=i + 1 + random.uniform(-.1, .2),
                               y=d + random.uniform(-.2, .2),
                               alpha=.3)
    else:
        pyplot.boxplot(vals, 0, 'rs', 1)
        for d in vals:
            pyplot.scatter(x=1 + random.uniform(-.1, .1),
                        y=d + random.uniform(-.2, .2),
                        alpha=.3)
    if yline is not None:
        if isinstance(yline, tuple):
            print(yline)
            assert len(yline)== 3,\
                "yline must either be a single number or a 3-lenght tuple"
            pyplot.axhline(yline[0], color="green", linestyle="--"),
            pyplot.axhline(yline[1], color="green"),
            pyplot.axhline(yline[2], color="green", linestyle="--"),
        else:
            pyplot.axhline(yline, color="green"),
    pyplot.title(title)
    fig.savefig(outpath)


def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values.

    # borrowed from http://code.activestate.com/recipes/511478-finding-the-percentile-of-the-values/
    @parameter N - is a list of values.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    N.sort()
    if not N:
        return None
    k = (len(N) - 1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c - k)
    d1 = key(N[int(c)]) * (k - f)
    return d0 + d1


def populate_subgraph_from_source(g, root, node_list, counter, debug=False):
    # counter for dbugging
    if debug:
        print("populating recursion depth %i" % counter)
    # starting from node, examit its neighbors
    for neigh in root.neighbor_list:
        # find that node in the node_list by
        for node in node_list:
            # ...checking name and orientation same, and not already in graph
            if node.name == neigh.name and \
               node.reverse_complimented == neigh.reverse_complimented and \
               node.name not in g.nodes():
                # if found, add that node and the appropriate edges
                g.add_node(node.name, cov=node.cov, length=node.length)
                g.add_edge(root.name, node.name)
                # rinse and repeat
                populate_subgraph_from_source(
                    g=g,
                    root=node,
                    node_list=node_list,
                    counter=counter + 1)


def reverse_populate_subgraph_from_source(g, root, node_list, counter, debug=False):
    # counter for dbugging
    if debug:
        print("populating rev recursion depth %i" % counter)
    # print(root)
    # look through all nodes
    for node in node_list:
        # print(node.name)
        if node.name in g.nodes():
            continue
        # and each of thats node's neighbors
        for neigh in node.neighbor_list:
            # that list our root as its neighbor (with right orientation)
            if neigh.name == root.name and \
               neigh.reverse_complimented == root.reverse_complimented:
                g.add_node(node.name, cov=node.cov, length=node.length)
                # we want to build the directionality opposite what it is
                # currently, so we make the nodes going from the root to the node
                g.add_edge(root.name, node.name)
                # rinse and repeat
                reverse_populate_subgraph_from_source(
                    g=g,
                    root=node,
                    node_list=node_list,
                    counter=counter + 1)


def make_rRNAs_dict(gff_list, gff_list_partial):
    solid16, solid23, solid5 = find_rRNA_from_gffs(gff_list, partial=False)
    partial16, partial23, partial5 = find_rRNA_from_gffs(gff_list_partial, partial=True)
    partial16 = [x for x in partial16 if x not in solid16]
    partial23 = [x for x in partial23 if x not in solid23]
    partial5  = [x for x in partial5 if x not in solid5]
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
    return rrnas


def check_rrnas_dict(rrnas, logger=None):
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


def remove_duplicate_nested_lists(l):
    """ removes duplicates from nested integer lists
    """
    uniq = set(tuple(x) for x in l)
    L = [ list(x) for x in uniq ]
    return L


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

    # make a list of node objects, a adjacency matrix M, and a DiGRaph object G
    node_list, M, G = parse_fastg(f=args.assembly_graph)

    if PLOT:
        plot_adjacency_matrix(
            M, node_order=None, partitions=[], colors=[],
            outpath=os.path.join(output_root, "full_adjacency_matrix.pdf"))
    #  write out the adjacency matrix
    with open(os.path.join(args.output, "tab.txt"), "w") as o:
        for line in M:
            o.write("\t".join([str(x) for x in line]) + "\n")
    # alt nodes
    dic = {int(x.name): [int(y.name) for y in x.neighbor_list] for x in node_list}


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
        evalue=10,
        kingdom="bac")
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

    # dict  of {gene: {partial: [nodes]; solid: [nodes]}}has keys of gense, where the value are dict
    rrnas = make_rRNAs_dict(gff_list, gff_list_partial)
    logger.debug(rrnas)
    check_rrnas_dict(rrnas, logger=logger)

    # this holds the {data} of the nodes keyed by their name
    nodes_data = dict(G.nodes(data=True))

    # add temporary edges; later we will weight them and redo the directionality if needed
    for node in node_list:
        for neigh in node.neighbor_list:
            G.add_edge(neigh.name, node.name)

    # get the depths of the big contigs
    depths_of_all_nodes, ave_depth_all_node, qs_all = get_depth_of_big_nodes(
        G, threshold=0)
    depths_of_big_nodes, ave_depth_big_node, qs_big = get_depth_of_big_nodes(
        G, threshold=4000)
    if PLOT:
        make_silly_boxplot(
            vals=[depths_of_all_nodes, depths_of_big_nodes],
            names=["All nodes", "Nodes > 4kb"],
            title="Depths of Nodes in Assembly Graph",
            outpath=os.path.join(output_root, "average_node_depths.pdf")
        )

    ########################################################################
    # collapse nodes where partial loci neighbor full-length
    collapsed = []
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
                # solidlen = dict(G.nodes(solid))['length']
                if part in G.neighbors(solid):
                    # print(G.edges(solid))
                    for d in G.edges(part):
                        if d[1] != solid and d[1] not in G.neighbors(solid):
                            # make a bi-directional graph for now
                        #############################################3
                            """ can we remake the node_list object to reflect these changes?
                            """
                            G.add_edge(solid , d[1])
                            G.add_edge(d[1], solid)
                    # G.remove_node(part)
                        #############################################3
                    these_collapsed.append(part)

        collapsed.extend(these_collapsed)
    logger.info("marked %i nodes for collapsing:", len(collapsed))
    logger.info(collapsed)

    ########  Reduce this graph to all nodes within 20kb if a 16S region
    interior_nodes = []
    border_nodes = []
    for i in rrnas["16S"]["solid"]:
        interior, border = neighborhood_by_length(G, i, cutoff=1000, ignored_nodes=[])
        print("16S node %i (close) neighbors", i)
        print(interior)
        print("16S node %i (border) neighbors", i)
        print(border)
        interior_nodes.extend(interior)
        border_nodes.extend(border)
    if rrnas["16S"]["solid"] != rrnas["23S"]["solid"]:
        for i in rrnas["23S"]["solid"]:
            interior, border = neighborhood_by_length(
                G, i, cutoff=1000, ignored_nodes=rrnas["16S"]["solid"])
            print("23S node %i (close) neighbors", i)
            print(interior)
            print("23S node %i (border) neighbors", i)
            print(border)
            interior_nodes.extend(interior)
            border_nodes.extend(border)
    valid_nodes = [x for y in [interior_nodes, border_nodes] for x in y]
    print(len(set(valid_nodes)))
    bad_nodes = set(G.nodes).symmetric_difference(set(valid_nodes))
    print("removing %i of %i nodes that aren't near rDNA" %(len(bad_nodes), len(G.nodes)))
    for node in bad_nodes:
        G.remove_node(node)
    if PLOT and args.plot_graphs:
        plot_G(
            G,
            rrnas["5S"]["solid"],
            rrnas["16S"]["solid"],
            rrnas["23S"]["solid"],
            outpath=os.path.join(args.output, "test_post_reduction_G.pdf"),
            outpath2=os.path.join(args.output, "test_post_reduction_G_linegraph.pdf"),
        )

    #######   identify paths between the 16S and 23S, if its less than say 500bp
    # connector_paths = []
    # for node16 in rrnas["16S"]["solid"]:
    #     for node23 in rrnas["23S"]["solid"]:
    #         connector_paths.extend(
    #             nx.all_simple_paths(G, node16, node23, cutoff=500)
    #         )
    # print("number of connector paths: %i" %len(connector_paths))


    #####################################################################################
    """
    So at this point, we have a poor implementation of the directionality of
    the graph, which we will need in a minute for path finding.  What we
    will do here is make a brand new subgraph for both the 5' and 3' regions
    of the operpn.

    lets try to add some edges to our graph  First we need to identify our
    center. For simple instaces (which is all we handle for now), there are
    two scenarios:
    1) 1 clearly defined cluster, where the three genes are all on one contig.
    2) 1 clearly defined contig for each gene.

    In all cases, we set the center of the graph to be the end of the 16S
    gene.  We assume that they are both pointing in the same direction.
    If that is not the case, we have a problem.

    .         center/break
    .           *
    ========|16S| ==|23S|=|5S|=================
    <-----------  ----------------------->
    in both these cases, we take the same approach.  Look for paths from
    (reverse complimented) 16S and from 23S
    """

    subset_node_list = []
    for n in node_list:
        if n.name in G.nodes():
            subset_node_list.append(n)
    subgraphs = {}

    for region in ["16S", "23S"]:
        REV = False
        # lets check the gff to see if the 16S is a positive or a negative hit:
        for line in gff_list:
            if region in line[8]:
                print(line)
                if line[6] == "-":
                    REV = True
                break
        # whos on first?
        # all this does is to make sure we take the oposite path from
        # the 16S nodes
        reverse_compliment_now = REV if region == "16S" else not REV
        node_region_data = extract_node_len_cov_rc(line[0].split(":")[0])

        # retrieve the starting node (16S or 23S) from the list of nodes,
        # making sure we get the one with the correct orientation
        init_node = None
        for N in subset_node_list:
            if N.name == int(node_region_data[0]) and N.reverse_complimented ==  REV:
                init_node = N
        if init_node is None:
            raise ValueError("node note found to initiate recursive subtree construction")
        logger.debug("Initial node")
        logger.debug(init_node)

        # add this initial node to a brand new DiGraph
        g = nx.DiGraph()
        g.add_node(int(node_region_data[0]), cov=int(node_region_data[2]), length=int(node_region_data[1]))

        # here is the path tracing algorithm.
        # We get all paths leading to the 16S, or all the paths away from 23S
        if region == "16S":
            populate_subgraph_from_source(g=g, root=init_node, node_list=subset_node_list, counter=1)
        else:
            reverse_populate_subgraph_from_source(g=g, root=init_node, node_list=subset_node_list, counter=1)

        logger.debug("nodes in %s subgraph", region)
        logger.debug(g.nodes())


        if PLOT and args.plot_graphs:
            logger.info("plotting reconstructed tree from %s node  %s" %\
                        (region, init_node.name))
            plot_G(
                g,
                rrnas["5S"]["solid"],
                rrnas["16S"]["solid"],
                rrnas["23S"]["solid"],
                outpath=os.path.join(args.output, "tiny%s.pdf" % region),
                outpath2=os.path.join(
                    args.output, "line_graph_tiny%s.pdf" %region),
            )
        if PLOT:
            plot_adjacency_matrix(
                make_adjacency_matrix(nx.convert.to_dict_of_dicts(g)), node_order=None,
                partitions=[], colors=[],
                outpath=os.path.join(
                    args.output,
                    "%s_subgraph_adjacency_matrix.pdf" % region))

        # count the paths going out from the 16S/23S
        out_paths_region_raw = []
        out_paths_region_sans_collapsed = []
        out_paths_region = []

        logger.debug("tips will have an out-degree of 0 on this reduced graph")
        logger.debug([g.out_degree(node) for node in g.nodes()])
        tips = [node for node in g.nodes() if g.out_degree(node) == 0]
        logger.debug("tips of subgraph:")
        logger.debug(tips)
        for noderegion in rrnas[region]["solid"]:
            for tip in tips:
                out_paths_region_raw.extend(nx.all_simple_paths(g, noderegion, tip))
        logger.info("number of raw paths to %s: %i" %(region, len(out_paths_region_raw)))

        """
        remember the collaping stuff we figured out before?  Lets remove
        those collapsable nodes from the paths now, and get rid of
        redundant paths.  We also take this opportunity to filter out
        paths that include tips not at the tip.  It happens, somehow...
        """
        for path in out_paths_region_raw:
            new_path = []
            internal_tip = False
            for i in path:
                if i not in collapsed:
                    new_path.append(i)
            # here we only add it the path doesnt contain tips (excluding the last item)
            for tip in tips:
                if tip in path[0: len(path) - 1]:
                    internal_tip = True
                else:
                    pass
            if not internal_tip:
                out_paths_region_sans_collapsed.append(new_path)

        # filter out duplicated paths:
        out_paths_region_sans_collapsed = remove_duplicate_nested_lists(
            out_paths_region_sans_collapsed)
        # now we have filtered, and some paths might not be long enough. Others, if the graph is cyclical, will be too long.  here, we trim and filter!
        for path in out_paths_region_sans_collapsed:
            filtered_path = [path[0]]
            filtered_length = 0
            internal_tip = False
            # skip first node, our root (either 16S or 23S)
            for i, node in enumerate(path):
                if i > 0:
                    node_length = nodes_data[node]["length"]
                    filtered_length = filtered_length + node_length
                    filtered_path.append(node)
                    if filtered_length > 1000:
                        break
            # if path does indeed pass the threshold
            if filtered_length > 1000:
                out_paths_region.append(filtered_path)

        # filter out duplicated paths, again:
        out_paths_region = remove_duplicate_nested_lists(out_paths_region)

        logger.debug(
            "Removed %i paths that contained collapsable nodes, etc" % \
            (len(out_paths_region_raw) - len(out_paths_region)))
        logger.info("%s Paths:\n", region)
        for i in out_paths_region:
            logger.info(i)

        all_region_path_nodes = [x for y in out_paths_region for x in y]
        set_region_path_nodes_normalized_depth = {}
        for i in set(all_region_path_nodes):
            set_region_path_nodes_normalized_depth[i] = nodes_data[i]['cov'] / all_region_path_nodes.count(i)

        if PLOT:
            make_silly_boxplot(
                vals=[x for x in set_region_path_nodes_normalized_depth.values()],
                outpath=os.path.join(output_root, "normalized_depths_of_%s_nodes.pdf" % region),
                title= "Normalized depths of nodes connected to %s" % region,
                yline=(qs_big[0], ave_depth_big_node, qs_big[2])
            )

        logger.debug("determining the depths of the %s paths" % region)
        all_region_path_depths = []
        for i, path in enumerate(out_paths_region):
            sublist = []
            for node in path:
                sublist.append(set_region_path_nodes_normalized_depth[node])
            all_region_path_depths.append(sublist)
        if PLOT:
            make_silly_boxplot(
                vals=all_region_path_depths,
                outpath=os.path.join(output_root, "normalized_depths_of_%s_exiting_paths.pdf" % region),
                title= "Depths of paths connected to %s" % region,
                yline=(qs_big[0], ave_depth_big_node, qs_big[2])
            )
        subgraphs[region] = {
            "raw_paths": out_paths_region_raw,
            "filtered_paths": out_paths_region,
            "graph": g
        }


    # interpret results
    n_upstream = len(subgraphs["16S"]["filtered_paths"])
    n_downstream = len(subgraphs["23S"]["filtered_paths"])
    logger.info("Paths leading to rDNA operon: %i", n_upstream)
    logger.info("Paths exiting  rDNA operon: %i", n_downstream)
    if n_upstream == n_downstream:
        logger.info("This indicates that there are at least %i rDNAs present",
                    n_upstream)
    else:
        logger.info("inconclusive results")
