#!/usr/bin/env python
"""
version 0.8.3
Minor version changes:
 - offloaded functions to utils
 - cleaned up old code gunk
 - applied WD-40 to the bad parts
 - added width option to include within regions
 - changed all mentions of width to within, cause width is not accurate
 - replace argument added
#TODO:
- set up logging
- Make this less awful
- Make this work with multiple scaffolds.  Should just mean
    wrapping certain looks in another for rec in record bit
- Maybe use percentage based exclusion, rather than a "within" arg?
Input:
- genbank file
- dictionary
- specific features : 16S, 5S
- upstream, downstream widths

Output:
-dir containing DNA fastas in their

"""
import re
import os
import pprint
#import random
#import sys
import csv
import subprocess
import datetime
import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from pyutilsnrw import utils3_5
from pyutilsnrw.utils3_5 import check_single_scaffold, get_genbank_seq, \
    get_genbank_record

pp = pprint.PrettyPrinter(indent=4)

def get_args():
    parser = argparse.ArgumentParser(description="This is used to extract regions " +
                                     " of interest based on supplied locus tags.")
    parser.add_argument("genbank_genome", help="Genbank file (WITH SEQUENCE)")
    parser.add_argument("clustered_loci", help="output from riboSelect")
    parser.add_argument("-f", "--feature", help="Feature, such as CDS,tRNA, " +
                        "rRNA; default: %(default)s",
                        default='rRNA', dest="feature",
                        action="store", type=str)
    parser.add_argument("-s", "--specific_features", help="colon:separated " +
                        "-- specific features to be grepped from product, " +
                        "such as 16S or tRNA-Ala; default: %(default)s",
                        default='16S:23S', type=str, dest="specific_features",
                        action="store")
    parser.add_argument("-w", "--within_feature_length",
                        help="bp's to include within the region; " +
                        "default: %(default)s",
                        default=0, dest="within", action="store",type=int)
    parser.add_argument("-m", "--minimum_feature_length",
                        help="if --replace, and a sequence is shorter than " +
                        " 2x --within_feature_length, --within will be " +
                        " modified so that only -m bp of sequnece are" +
                        "turned to N's " +
                        "default: %(default)s",
                        default=100, dest="minimum", action="store",type=int)
    parser.add_argument("-l", "--flanking_length",
                        help="length of flanking regions, can be colon-" +
                        "separated to give separate upstream and " +
                        "downstream flanking regions; default: %(default)s",
                        default='700',type=str, dest="flanking")
    parser.add_argument("-r", "--replace",
                        help="replace sequence with N's; default: %(default)s",
                        default=False, action="store_true", dest="replace")
    parser.add_argument("-v", "--verbose",
                        help="verbose output (its ugly) default: %(default)s",
                        default=False, action="store_true", dest="verbose")
    parser.add_argument("-o", "--output",
                        help="output directory; default: %(default)s",
                        default=os.getcwd(),
                        type=str, dest="output")
    args = parser.parse_args()
    return(args)


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
            if (feat.qualifiers.get("locus_tag")[0] in locus_tag_list):  # and\
               # (feat.type == args.feature):
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



if __name__ == "__main__":
    args = get_args()
    # specific_features = args.specific_features.split(":")
    # feature_regs = args.feature_regions.split(":")
    # print(args)
    print("Usage:\n{0}\n".format(str(" ".join([x for x in sys.argv]))))

    # try:
    #     feature_regs_len = [int(x) for x in args.feature_regions_lengths.split(":")]
    # except ValueError:
    #     print("feature_regions_lengths must be an integer, or integers separated" +
    #           "by colons, such as 400:33")
    #     sys.exit(1)
    # if not len(specific_features) == len(feature_regs) == len(feature_regs_len):
    #     print("features, feature_refions, and feature_region_lengths must be equal")
    #     sys.exit(1)
    date = str(datetime.datetime.now().strftime('%Y%m%d'))
    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    clusteredDict = parse_clustered_loci_file(args.clustered_loci)
    genome_sequence = get_genbank_seq(args.genbank_genome)
    genome_record = get_genbank_record(args.genbank_genome)
    regions = []
    for i in clusteredDict.keys():
        locus_tag_list = clusteredDict[i]
        coord_list = extract_coords_from_locus(genome_record,
                                               locus_tag_list=locus_tag_list,
                                               verbose=args.verbose)
        regions.append(stitch_together_target_regions(genome_sequence,
                                                      coords=coord_list,
                                                      within=args.within,
                                                      minimum=args.minimum,
                                                      flanking=args.flanking,
                                                      replace=args.replace,
                                                      verbose=args.verbose))
    output_index = 1
    for i in regions:
        filename = str("region_" + str(output_index))
        i = str(i)
        # print(i)
        # print(len(i))
        with open(os.path.join(args.output,
                               str(date + "_" + filename + "_riboSnag.fasta")),
                  "w") as outfile:
            SeqIO.write(SeqRecord(Seq(i,IUPAC.IUPACAmbiguousDNA()),
                                  id = str("{0}_riboSnag_{1}_flanking_{2}_within".format(
                                      output_index, args.flanking, args.within)),
                                  description=""), outfile, "fasta")
            outfile.write('\n')
            output_index = output_index + 1
