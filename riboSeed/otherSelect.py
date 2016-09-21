#!/usr/bin/env python
"""
With vesion 0.8.0 of riboseed, we moved away from the flexibility of getting any region from a genome.  This is a stand in to create a list of other loci of interest. 

version 0.8.1
Minor version changes:
 - This is the first instance of otherSelect.py and riboSelect.py.
#TODO:
- Make this less awful
- No seriously, this is awful
Input:
- genbank file
- location: cds, gene, or rRNA
- specific features : 16S, 5S
- specific regions : upstream, downstream, etc
- region lengths: int

Output:
-DNA fasta

USAGE:
 $ python riboSnag.py genome.gbk --features rRNA --specific_features 16S:23S --feature_regs upstream:downstream --feature_regs_len 200:400
     --within_feature_length 50  --output /path/to/output
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
pp = pprint.PrettyPrinter(indent=4)

#%%


def get_args(DEBUG=False):
    # dummy_args = argparse.Namespace(input_genome_path="NC_011751.1.gb",
    #                                 feature="rRNA",
    #                                 pre_specific_features="16S:23S",
    #                                 pre_feature_regs="upstream:downstream",
    #                                 pre_feature_regs_len="200:400",
    #                                 within=50,
    #                                 output=os.path.expanduser("~/Desktop/"))
    # if DEBUG:  # ie, if running interactively
    #     print("Using Dummy Arguments")
    #     return(dummy_args)
    parser = argparse.ArgumentParser(description="This is used to extract regions of interest based\
                                                  on coordinates and feature types.")
    parser.add_argument("genbank_genome", help="Genbank file (WITH SEQUENCE)")
    parser.add_argument("-f", "--feature", help="Feature, such as CDS,tRNA, rRNA; default: %(default)s",
                        default='rRNA', type=str)
    parser.add_argument("-s", "--specific_features", help="colon:separated -- specific features\
                         to be grepped from product, such as 16S or tRNA-Ala; default: %(default)s",
                        default='16S:23S', type=str)
    parser.add_argument("-r", "--feature_regions", help="colon:separated -- upstream, \
                         upstream_within, downstream, downstream_within, identity, identityN; default: %(default)s",
                        default="upstream_within:downstream_within", type=str)
    parser.add_argument("-l", "--feature_regions_lengths", help="colon:separated -- length of\
                         region; default: %(default)s",  default="100:100", type=str)
    parser.add_argument("-w", "--within_feature_length", help="if using up/downstream_within, \
                         where you return some of the feature itself, this dictates how many bp; default: %(default)s",
                        default=50, type=int)
    parser.add_argument("-o", "--output", help="output directory; default: %(default)s", default=os.getcwd(),
                        type=str, dest="output")
    # parser.add_argument("--DEBUG", help="run on test data; default: %(default)s", default=False,
    #                     action="store_true", dest="DEBUG")
    args = parser.parse_args()
    # if args.DEBUG:  # ie, if running debug from cline
    #     print("Using Dummy Arguments")
    #     return(dummy_args)
    # else:
    return(args)


def check_single_scaffold(input_genome_path):
    """Test for single scaffold
    """
    print("testing for multiple scaffolds...")
    counter = -1  # if all goes well, returns 0, else returns 1 or more or -1
    for line in open(input_genome_path, "r"):
        if re.search("ORIGIN", line) is not None:
            counter = counter + 1
    return(counter)


def get_genbank_record(input_genome_path):
    """reads the FIRST record only from a genbank file; will probably only work for first scaffold
    """
    print("Reading genbank file...")
    with open(input_genome_path) as input_genome_handle:
        genome_seq_record = next(SeqIO.parse(input_genome_handle, "genbank"))
    if genome_sequence[0: 100] == str("N" * 100):
        print("Careful: the first 100 nucleotides are N's; did you download the full .gb file?")
    return(genome_seq_record)


def get_genbank_seq(input_genome_path):
    """get the sequence from the FIRST record only in a genbank file
    """
    print("fetching nucleotide sequence from genbank file...")
    with open(input_genome_path) as input_genome_handle:
        genome_seq_record = next(SeqIO.parse(input_genome_handle, "genbank"))
    return(genome_seq_record.seq)


# def extract_coords(genome_seq_record, feature="rRNA", verbose=True):
#     """the magic (part 1) of this whole thing
#     """
#     loc_number = 0  # index for hits
#     loc_list = []  # recipient structure
#     for feat in genome_seq_record.features:
#         if feat.type != feature:  # skip if not features of interest
#             pass
#         else:
#             #  SeqIO subtracts one from start coord to make it pythonic. the +1 below undoes that
#             coords = [feat.location.start.position + 1, feat.location.end.position]
#             strand = feat.strand
#             product = feat.qualifiers.get("product")
#             if "locus_tag" in feat.qualifiers:
#                 locus_tag = feat.qualifiers.get("locus_tag")
#             else:
#                 locus_tag = "none"
#             loc_list.append([loc_number, coords, strand, product, locus_tag])
#             loc_number = loc_number + 1
#     if not loc_number > 0:
#         raise ValueError("no hits found!")
#         sys.exit(1)
#     if verbose:
#         print("Here are the detected region, strand, product, locus tag, \
#                and subfeatures of the results:") 
#         pp.pprint(loc_list)
#     return(loc_list)


def extract_coords_from_locus(genome_seq_record, locus_tag_list=[], verbose=True):
    """given a list of locus_tags, return a list of loc_number,coords, strand, product
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
                locus_tag = feat.qualifiers.get("locus_tag")
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

def pure_python_kmeans(data, centers=3):
    """giveb 1d data and number of centers, returns a
    csv with the data and cluster, and LP's disapointment
    """
    with open(os.path.join(os.getcwd(), "list.csv"), "w") as f:
        for i in data:
            f.write("".join([str(i), "\n"]))
    
    with open(os.path.join(os.getcwd(), "km_script.R"), "w") as f:
        rcmds= ["centers<-{0}".format(centers),
                "data<-read.csv('list.csv', header=F, col.names=c('index'))", 
                "km <- kmeans(data[,1], nstart=100, iter.max=100, centers=centers)",
                "data[,2]<-km[1]",
                "data<- data[order(data[,1]),  ]",
                "write.csv(data, 'list.csv',  row.names=F)"]
        for i in rcmds:
                f.write(i)
                f.write('\n')
    subprocess.call("Rscript km_script.R", shell=sys.platform!='win32')

def get_filtered_locus_tag_dict(genome_seq_record, feature="rRNA", verbose=True):
    """returns dictionary of index:locus_tag id pairs for all 
    "feature"  entries.  This then gets clustered
    requires having locus tag in your genbank file.  Non-negitable
    """
    loc_number = 0 # counter
    locus_tag_dict = {}  # recipient structure
    for feat in genome_seq_record.features:
        try:
            locustag = feat.qualifiers.get("locus_tag")[0]
            locus_tag_dict[loc_number] = (locustag, feat.type) 
            loc_number = loc_number + 1
        except TypeError:
            pass
    if len(locus_tag_dict) < 1:
            raise ValueError("no locus tags found!")
    print("filtering by feature of interest")
    filtered = {k: v[0] for k, v in locus_tag_dict.items() if v[1] == feature}
    if verbose: 
        for key in sorted(filtered):
                print("%s: %s" % (key, filtered[key]))
    return(filtered)


def stitch_together_target_regions(genome_sequence, coords, flanking=400, replace=True):
    """ 
    given a list from get_coords, usually of length 3 (16,5,and 23 rRNAs),
    return a string with the sequence of the region, replacing coding sequences
    with N's (or not, replace=False), and including the flanking regions 
    upstream and down.
    """
    region = ''
    for i in range(0, len(coords)):
        if i == 0:  # if first in list, start with flanking prefix
            start = coords[i][1][0] - flanking
            end = coords[i][1][0]
            region = "".join([region, str(genome_sequence[start:end])])
        if replace:
            start = coords[i][1][0]
            end = coords[i][1][1]
            region = "".join([region, str("N" * (end - start + 1))])
        else:
            region = "".join([region, str(genome_sequence[coords[i][1][0]: coords[i][1][1]])])               
        if i == (len(coords)-1):  # if last in list, tag on flanking suffix 
            start = coords[i][1][1]
            end = coords[i][1][1] + flanking
            region ="".join([region,  str(genome_sequence[start:end])])
        else:  # if not, get intergenic region, if it exists
            start = coords[i][1][1] + 1
            end = coords[i + 1][1][0] - 1
            region ="".join([region,  str(genome_sequence[start:end])])
    return(region)

def get_region_near_loc(genome_sequence, loc_list, feature_list, feature_regions,
                        feature_regs_len, within):
    """ The magic, part 2
    upstream, downstream, identity, identityN for reiong repped as N's
    TODO make this fail if it tries to reference outside of bounds
    """
    feature_region_options = ["upstream", "upstream_within", "downstream",
                              "downstream_within", "identity",
                              "identityN"]
    if not len(feature_list) == len(feature_regions) == len(feature_regs_len):
        raise ValueError("you supplied  as a feature list; length of feature_regions must \
                          have same length to match up properly")
    for i in range(0, len(feature_list)):
        # check feature
        if not feature_regions[i] in feature_region_options:
            raise ValueError(str(
                "{0} not in availible feature regions below:\n{1}".format(
                    feature_regions[i], feature_region_options)))
        for j in loc_list:
            if re.search(feature_list[i], str(j[3])):
                fiveprime = None
                threeprime = None
                dna = None
                if feature_regions[i] == "upstream" and j[2] == 1:
                    fiveprime = j[1][0] - feature_regs_len[i] - 1  # -1 accounts for 0 index
                    threeprime = j[1][0]
                    dna = genome_sequence[fiveprime:threeprime]
                    desc = "_from_%i_to_%i" % (fiveprime + 1, threeprime)
                elif feature_regions[i] == "upstream" and j[2] == -1:
                    fiveprime = j[1][1] + 1
                    threeprime = j[1][1] + feature_regs_len[i]
                    dna = genome_sequence[fiveprime:threeprime].reverse_complement().lower()
                    desc = "_from_%i_to_%i, revcomplemented" % (fiveprime + 1, threeprime)
                elif feature_regions[i] == "upstream_within" and j[2] == 1:
                    fiveprime = j[1][0] - feature_regs_len[i] - 1
                    threeprime = j[1][0] + within
                    dna = genome_sequence[fiveprime:threeprime]
                    desc = "_from_%i_to_%i" % (fiveprime + 1, threeprime)
                elif feature_regions[i] == "upstream_within" and j[2] == -1:
                    fiveprime = j[1][1] + 1 - within
                    threeprime = j[1][1] + feature_regs_len[i]
                    dna = genome_sequence[fiveprime:threeprime].reverse_complement().lower()
                    desc = "_from_%i_to_%i, revcomplemented" % (fiveprime + 1, threeprime)
                elif feature_regions[i] == "downstream" and j[2] == 1:
                    fiveprime = j[1][1] + 1
                    threeprime = j[1][1] + feature_regs_len[i]
                    dna = genome_sequence[fiveprime:threeprime]
                    desc = "_from_%i_to_%i" % (fiveprime + 1, threeprime)
                elif feature_regions[i] == "downstream" and j[2] == -1:
                    fiveprime = j[1][0] - feature_regs_len[i] - 1
                    threeprime = j[1][0]
                    dna = genome_sequence[fiveprime:threeprime].reverse_complement().lower()
                    desc = "_from_%i_to_%i, revcomplemented" % (fiveprime + 1, threeprime)
                elif feature_regions[i] == "downstream_within" and j[2] == 1:
                    fiveprime = j[1][1] + 1 - within
                    threeprime = j[1][1] + feature_regs_len[i]
                    dna = genome_sequence[fiveprime:threeprime]
                    desc = "_from_%i_to_%i" % (fiveprime + 1, threeprime)
                elif feature_regions[i] == "downstream_within" and j[2] == -1:
                    fiveprime = j[1][0] - feature_regs_len[i] - 1
                    threeprime = j[1][0] + within
                    dna = genome_sequence[fiveprime:threeprime].reverse_complement().lower()
                    desc = "_from_%i_to_%i, revcomplemented" % (fiveprime + 1, threeprime)
                elif feature_regions[i] == "identity" and j[2] == 1:
                    fiveprime = j[1][0] - 1
                    threeprime = j[1][1]
                    dna = genome_sequence[fiveprime:threeprime]
                    desc = "_from_%i_to_%i" % (fiveprime + 1, threeprime)
                elif feature_regions[i] == "identity" and j[2] == -1:
                    fiveprime = j[1][0] - 1
                    threeprime = j[1][1]
                    dna = genome_sequence[fiveprime:threeprime].reverse_complement().lower()
                    desc = "_from_%i_to_%i, revcomplemented" % (fiveprime + 1, threeprime)
                elif feature_regions[i] == "identityN":
                    fiveprime = j[1][0] - 1
                    threeprime = j[1][1]
                    dna = str("N" + ("N" * (fiveprime - threeprime)))
                    desc = "_from_%i_to_%i" % (fiveprime + 1, threeprime)
                else:
                    raise ValueError("Why isnt %s caught earlier" %
                                     feature_regions[i])
                j.append(dna)
                j.append(desc)
            else:
                continue
    rm_list = []  # remove entries without sequences
    for i in loc_list:
        try:
            len(i[5])
        except IndexError:
            rm_list.append(i)
    trimmed_list = [x for x in loc_list if x not in rm_list]
    return(trimmed_list)



#def merge_neighboring_sequences(coords_with_seq, window_width=500):
#    group_number = 10
#    coords_with_seq[0].append(10)
#    for i in range(1, len(coords_with_seq)):
#        if coords_with_seq[i][1][0] - coords_with_seq[i-1][1][1] < window_width:
#            coords_with_seq[i-1].append(group_number)
#            print(group_number)
#        else:
#            group_number = group_number + 10
#        ## take care of last one; below statement should always be true in a perfect world
#    if coords_with_seq[len(coords_with_seq)-1][1][0] - coords_with_seq[len(coords_with_seq)-2][1][1] < window_width :
#        coords_with_seq[len(coords_with_seq)-1].append(group_number)
#    for i in range(10, group_number+10, 10):
#        coords_gh = []
#        ncoords = 0
#        for j in coords_with_seq:
#            if j[-1] == i:
#                coords_gh.append(j[1])
#                ncoords = ncoords + 1
#        print(coords_gh)
#        print("hi")
#merge_neighboring_sequences(coords_with_seq, windxoow_width=500)
#%%

if __name__ == "__main__":
    raise("this doesnt work yet")
# make lists if colon:delimeted
    args = get_args(DEBUG=False)
    specific_features = args.specific_features.split(":")
    feature_regs = args.feature_regions.split(":")
    print(args)
    try:
        feature_regs_len = [int(x) for x in args.feature_regions_lengths.split(":")]
    except ValueError:
        print("feature_regions_lengths must be an integer, or integers separated" +
              "by colons, such as 400:33")
        sys.exit(1)
    if not len(specific_features) == len(feature_regs) == len(feature_regs_len):
        print("features, feature_refions, and feature_region_lengths must be equal")
        sys.exit(1)
    date = str(datetime.datetime.now().strftime('%Y%m%d'))
    if check_single_scaffold(args.genbank_genome):
        raise ValueError("You really should only give this genbank files with one \
                         scaffold for now. Problem with that? Fork off.")
    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    genome_sequence = get_genbank_seq(args.genbank_genome)
    genome_record = get_genbank_record(args.genbank_genome)
    # TODO ofload this part til the #******# to  riboSelect.py
    lociDict = get_filtered_locus_tag_dict(genome_record)
    # for each specific feature, get number of occurances
    nfeatures_occur = []
    for i in specific_features:
        unique=[]
        for k,v in lociDict.items():
            t = v.find(i)
            unique.append(t)
        nfeatures_occur.append(len([x for x in unique if x!=-1]))
    print(str("number of occuraces of each specific feature: {0}; using " +
              " {1} clusters").format(nfeatures_occur, min(nfeatures_occur) ))
    pure_python_kmeans(lociDict.keys(), centers=min(nfeatures_occur))
    ## csv to dict
    with open('list.csv', mode='r') as infile:
        reader = csv.reader(infile)
        next(reader, None)  # skip the headers
        indexClusterDict = dict((rows[0],rows[1]) for rows in reader)
    ##
    print(lociDict)
    print(indexClusterDict)
    clusteredDict={}
    for k,v  in indexClusterDict.items():
        clusteredDict.setdefault(v,[]).append(lociDict[int(k)])
    print(clusteredDict)
    #******#
    print(clusteredDict["1"])
    regions = []
    for i in clusteredDict.keys():
        locus_tag_list = clusteredDict[i]
        coord_list = extract_coords_from_locus(genome_record,
                                               locus_tag_list=locus_tag_list,
                                               verbose=True)
        regions.append(stitch_together_target_regions(genome_sequence,
                                                     coords=coord_list,
                                                     flanking=400,
                                                     replace=True))
    output_index = 1
    for i in regions:
        filename = str("region_" + str(output_index))
        with open(os.path.join(args.output,
                               str(date + "_" + filename + "_riboSnag.fasta")),
                  "w") as outfile:
            SeqIO.write(SeqRecord(Seq(i,IUPAC.IUPACAmbiguousDNA()),
                                  id = str(str(output_index) + "_riboSnag"),
                                  description=""), outfile, "fasta")
            outfile.write('\n')
            output_index = output_index + 1
