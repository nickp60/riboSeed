#!/usr/bin/env python
"""
version 0.8.5
Minor version changes:
 - super-prelim multiple cscaffold handling
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
import csv
import subprocess
import datetime
import time
import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

#from pyutilsnrw import utils3_5
from pyutilsnrw.utils3_5 import get_genbank_seq, get_genbank_record, \
    set_up_logging


def get_args():
    parser = argparse.ArgumentParser(description="Use to extract regions " +
                                     "of interest based on supplied locus " +
                                     " tags.")
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
                        default=0, dest="within", action="store", type=int)
    parser.add_argument("-m", "--minimum_feature_length",
                        help="if --replace, and a sequence is shorter than " +
                        " 2x --within_feature_length, --within will be " +
                        " modified so that only -m bp of sequnece are" +
                        "turned to N's " +
                        "default: %(default)s",
                        default=100, dest="minimum", action="store", type=int)
    parser.add_argument("-l", "--flanking_length",
                        help="length of flanking regions, can be colon-" +
                        "separated to give separate upstream and " +
                        "downstream flanking regions; default: %(default)s",
                        default='700', type=str, dest="flanking")
    parser.add_argument("-r", "--replace",
                        help="replace sequence with N's; default: %(default)s",
                        default=False, action="store_true", dest="replace")
    parser.add_argument("-p", "--per_contig",
                        help="if genome is not small or nearly completed, " +
                        "use --per_contigs. This searches for loci cluster " +
                        "on a per-contig  basis as opposed to globally " +
                        "searching for loci; default: %(default)s",
                        default=False, action="store_true", dest="per_contig")
    parser.add_argument("-v", "--verbosity", dest='verbosity', action="store",
                        default=2, type=int,
                        help="1 = debug(), 2 = info(), 3 = warning(), " +
                        "4 = error() and 5 = critical(); default: %(default)s")
    parser.add_argument("-o", "--output",
                        help="output directory; default: %(default)s",
                        default=os.getcwd(),
                        type=str, dest="output")
    parser.add_argument("--clobber",
                        help="overwrite previous output files" +\
                        "default: %(default)s", action='store_true',
                        default=False, dest="clobber")
    args = parser.parse_args()
    return(args)


def parse_clustered_loci_file(file, logger=None):
    """changed this to return a list, cause multiple clusters
    can be on the smae contig, which makes for duplicated
    dict keys, and those dont work
    """
    if logger is None:
        raise ValueError("logging must be used!")
    if not os.path.exists(file):
        logger.error("Cluster File not found!")
        sys.exit(1)
    clusters = []
    try:
        with open(file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    print("look, a coment!")
                seqname = line.strip("\n").split(" ")[0]
                clusters.append([seqname,
                                 [line.strip("\n").split(" ")[1].split(":")]])
    except:
        logger.error("Cluster file could not be parsed!")
        sys.exit(1)
    if len(clusters) == 0:
        logger.error("Cluster file could not be parsed!")
        sys.exit(1)
    return(clusters)


def extract_coords_from_locus(genome_seq_records, locus_tag_list=[],
                              feature="rRNA", verbose=True, logger=None):
    """given a list of locus_tags, return a list of
    loc_number,coords, strand, product
    """
    if logger is None:
        raise ValueError("logging must be used!")
    loc_number = 0  # index for hits
    loc_list = []  # recipient structure
    for record in genome_seq_records:
        logger.debug("searching {0} for loci in this list: {1}".format(
            record.id, locus_tag_list))
        for feat in record.features:
            if not feat.type in feature:
                continue
            logger.debug("found {0} in the following feature : \n{1}".format(
                feature, feat))

            try:
                locus_tag = feat.qualifiers.get("locus_tag")[0]
            except:
                logger.error(str("found a feature ({0}), but there is no" +\
                                 "locus tag associated with it!").format(
                                     feat))
                sys.exit(1)

            if locus_tag in locus_tag_list:
                #  SeqIO makes coords 0-based; the +1 below undoes that
                coords = [feat.location.start.position + 1,
                          feat.location.end.position]
                strand = feat.strand
                product = feat.qualifiers.get("product")
                # locus_tag = feat.qualifiers.get("locus_tag")[0]
                loc_list.append([loc_number, coords, strand,
                                 product, locus_tag, record.id])
                loc_number = loc_number + 1
            else:
                pass
    if not loc_number > 0:
        logger.error("no hits found in any record! Double " +
                     "check your genbank file")
        sys.exit(1)
    logger.debug("Here are the detected region,coords, strand, product, " +
                "locus tag, subfeatures aand sequence id of the results:")
    logger.info(loc_list)
    return(loc_list)



def get_genbank_seq_matching_id(recordID, genbank_record_list):
    """ given a list of loci and genbank records, return sequence of
    genbank record that has all the loci.
    If on different sequences, return error
    """
    if logger is None:
        raise ValueError("logging must be used!")
    for record in genbank_record_list:
        if recordID == record.id:
            return(record.seq)
        else:
            pass
    print("no seqeuence found matching record id!")
    sys.exit(1)


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
    if logger is None:
        logger.error("Must have logger for this function")
        sys.exit(1)
    try:
        flank = [int(x) for x in flanking.split(":")]
        if len(flank) == 1:  # if only one value use for both up and downstream
            flank.append(flank[0])
        assert(len(flank) == 2)
    except:
        logger.error("Error parsing flanking value; must either be " +
                     " integer or two colon-seapred integers")
        sys.exit(1)
    #TODO : make this safer
    smallest_feature = min([y[1] - y[0] for y in [x[1] for x in coords]])
    if smallest_feature < (minimum):
        raise ValueError("invalid minimum! cannot exceed half of smallest " +
                         "feature, which is {0} in this case".format(
                             smallest_feature))
    if verbose:
        for i in coords:
            logger.info(str(i))
    #  This works as long as coords are never in reverse order
    # TODO check coords are increasing
    global_start = min([y[0] for y in [x[1] for x in coords]]) - flank[0]
    # if start is negative, just use 0, the beginning of the sequence
    if global_start < 0:
        logger.warning("Caution! Cannot retrieve full flanking region, as " +\
                       "the 5' flanking region extends past start of sequence")
        global_start = 0
    global_end = max([y[1] for y in [x[1] for x in coords]]) + flank[1]
    if global_end > len(genome_sequence):
        logger.warning("Caution! Cannot retrieve full flanking region, as " +\
                       "the 3' flanking region extends past end of sequence")
        global_end = len(genome_sequence)
    full_seq = genome_sequence[global_start : global_end]
    seq_with_ns = str(full_seq)
    #
    # loop to mask actual coding regions with N's
    #
    if replace:
        for i in coords:
            region_length = i[1][1] - i[1][0] - (2 * within)
            # if dealing with short sequences
            if (i[1][1] - i[1][0]) < (2 * within):
                # set within to retain minimum sequence length
                this_within = int((i[1][1] - i[1][0] - minimum) / 2)
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
                              str("N" * region_length) +
                              seq_with_ns[rel_end: ])

    #
    #
    try:
        # make sure the sequence is correct length, corrected for zero-index
        assert(global_end - global_start, len(full_seq))
        # make sure replacement didnt change seq length
        assert(len(full_seq), len(seq_with_ns))
    except:
        logger.error("There appears to be an error with how the seqeuence " +
                     "coordinates are being calculated!")
    logger.info(str("\nexp length {0} \nact length {1}".format(
        global_end - global_start, len(full_seq))))
    if verbose:
        lb = 60
        for i in range(0, int(len(seq_with_ns) / lb)):
            print(str(full_seq[i * lb: lb + (i * lb)] + "\n"))
            print(str(seq_with_ns[i * lb: lb + (i * lb)] + "\n"))
            print("\n")
    seqrec = SeqRecord(Seq(seq_with_ns, IUPAC.IUPACAmbiguousDNA()),
                       id=str(coords[0][5] + "_" +
                              str(global_start) +
                              ".." + str(global_end)))
    return(seqrec)


if __name__ == "__main__":
    args = get_args()
    output_root = os.path.abspath(os.path.expanduser(args.output))
    # Create output directory only if it does not exist
    try:
        os.makedirs(args.output)
    except FileExistsError:
        print("#Selected output directory %s exists" %
              args.output)
        if not args.clobber:
            print("exiting")
            sys.exit(1)
        else:
            print("# continuing, and risking potential loss of data")
    logger = set_up_logging(verbosity=args.verbosity,
                            outfile=str("%s_riboSnag_log.txt" %
                                        os.path.join(output_root,
                                                     time.strftime("%Y%m%d%H%M"))),
                            name=__name__)

    print("Usage:\n{0}\n".format(str(" ".join([x for x in sys.argv]))))
    date = str(datetime.datetime.now().strftime('%Y%m%d'))
    # if not os.path.isdir(args.output):
    #     os.mkdir(args.output)
    clusteredList = parse_clustered_loci_file(args.clustered_loci,
                                              logger=logger)
    # genome_sequences = get_genbank_seq(args.genbank_genome, first_only=False)
    genome_records = get_genbank_record(args.genbank_genome, first_only=False,
                                        logger=logger)
    # print(genome_sequences)
    regions = []
    logger.info("clustered loci list:")
    logger.info(clusteredList)
    for i in clusteredList:
        locus_tag_list = i[1][0]
        recID = i[0]
        genbank_sequence = \
            get_genbank_seq_matching_id(recordID=recID,
                                        genbank_record_list=genome_records)
        coord_list = extract_coords_from_locus(genome_records,
                                               locus_tag_list=locus_tag_list,
                                               verbose=True, logger=logger)
        regions.append(stitch_together_target_regions(genbank_sequence,
                                                      coords=coord_list,
                                                      within=args.within,
                                                      minimum=args.minimum,
                                                      flanking=args.flanking,
                                                      replace=args.replace,
                                                      verbose=False,
                                                      logger=logger))
    logger.debug(regions)
    output_index = 1
    for i in regions:
        filename = str("region_" + str(output_index))
        with open(os.path.join(args.output,
                               str(date + "_" + filename + "_riboSnag.fasta")),
                  "w") as outfile:
            # i.description = str("{0}_riboSnag_{1}_flanking_{2}_within".format(
            #                           output_index, args.flanking, args.within))
            SeqIO.write(i, outfile, "fasta")
            outfile.write('\n')
            output_index = output_index + 1
