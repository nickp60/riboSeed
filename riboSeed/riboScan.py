#!/usr/bin/env python3
#-*- coding: utf-8 -*-

"""
"""

import argparse
import sys
import time
import glob
import re
import shutil
import subprocess
import logging
import os
import traceback
# import gffutils

from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.Alphabet import IUPAC

# need this line for unittesting
sys.path.append(os.path.join('..', 'riboSeed'))
from pyutilsnrw.utils3_5 import set_up_logging

# --------------------------- classes --------------------------- #


# --------------------------- methods --------------------------- #


def get_args():  # pragma: no cover
    """#TODO:     for cli mods:
    http://stackoverflow.com/questions/18025646/
         python-argparse-conditional-requirements
    make this able to handle different library types such as two unpaired runs
    """
    parser = argparse.ArgumentParser(
        description="Given de novo and de fere novo contigs files, a " +
        "misjoined de fere novo contig name,  and a colon:separated " +
        "list of de novo contig names, replace the offending contig with " +
        "the de novo contig(s) ",
        add_help=False)  # to allow for custom help
    parser.add_argument("contigs_dir", action="store",
                        help="multifasta containing de fere novo contigs")

    parser.add_argument("ext", action="store",
                        help="name of the bad contig")
    # taking a hint from http://stackoverflow.com/questions/24180527
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)

    # had to make this faux "optional" parse so that the named required ones
    # above get listed first
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-k", "--kingdom", dest='kingdom',
                          action="store",
                          choices=["bac", "euk"],
                          help="whether to look for eukaryotic or "+
                          "bacterial rDNA; " +
                          "default: %(default)s", default="bac",
                          type=str)
    optional.add_argument("-t", "--id_thresh", dest='id_thresh',
                          action="store", type=float,
                          help="; " +
                          "default: %(default)s", default=0.5,
                          )
    optional.add_argument("-b", "--barrnap_exe", dest='barrnap_exe',
                          action="store",
                          help="; " +
                          "default: %(default)s", default="barrnap",
                          type=str)
    optional.add_argument("-s", "--seqret_exe", dest='seqret_exe',
                          action="store",
                          help="; " +
                          "default: %(default)s", default="seqret",
                          type=str)
    optional.add_argument("-v", "--verbosity", dest='verbosity',
                          action="store",
                          default=2, type=int, choices=[1, 2, 3, 4, 5],
                          help="Logger writes debug to file in output dir; " +
                          "this sets verbosity level sent to stderr. " +
                          " 1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
                          "default: %(default)s")
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    args = parser.parse_args()
    return args


def last_exception():
    """ Returns last exception as a string, or use in logging.
    stolen verbatim from pyani
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


def parse_fasta_header(first_line):
    if not first_line.startswith(">"):
        raise ValueError("A valid fasta must start with a '>'!")
    if first_line.startswith(">gi"):
        accession = re.search(r">gi\|.*?\|.*?\|(.*?)\|",
                              first_line).groups()[0]
    else:
        accession = re.search(r">(.*?) .*", first_line).groups()[0]
    if accession is None:
        raise ValueError(
            "unable to extract accession from first line in file: \n%s".format(
                first_line))
    return accession


def make_barrnap_cmd(infasta, outgff, exe, thresh, kingdom):
    assert shutil.which(exe) is not None, "barrnap executable not found!"
    assert thresh > 0 and thresh < 1, "Thresh must be between 0 and 1!"
    cmd = "{0} -kingdom {1} {2} --reject {3} > {4}".format(
        shutil.which(exe),
        kingdom,
        infasta,
        thresh,
        outgff)
    return cmd


def add_locus_tags_to_gff(gff, acc):
    LOCUS = 0
    gff_list = []
    with open(gff, 'r') as g:
        for idx, line in enumerate(g):
            if idx == 0:
                gff_list.append(line)
            else:
                gff_list.append(
                    re.sub("Name=",
                           "locus_tag={0}_{1};Name=".format(acc, LOCUS),
                           line))
                LOCUS = LOCUS + 1

    new_gff = str(os.path.splitext(gff)[0] + "_tagged.gff")
    with open(new_gff, 'w') as outgff:
        for l in gff_list:
            print(str(l + "\n"))
            outgff.write(str(l + "\n"))
    return(new_gff)


def combine_gbs(finalgb, gb_list):
    with open(finalgb, 'w') as outfile:
        for fname in gb_list:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)


def append_accession_and_version(accession, ingb, finalgb):
    with open(finalgb, 'w') as outfile:
        with open(ingb) as infile:
            for idx, line in enumerate(infile):
                if idx == 1:
                    repline = "ACCESSION   {0}\nVERSION     {0}\n{1}".format(
                        accession, line)
                    outfile.write(repline)
                else:
                    outfile.write(line)


def make_seqret_cmd(exe, outgb, ingff, infasta):
    assert shutil.which(exe) is not None, "seqret executable not found!"
    cmd = str(
        "{0} -sequence {1} -feature -fformat gff3 -fopenfile {2} " +
        "-osformat genbank -auto  -outseq {3}"
    ).format(
        shutil.which(exe),
        infasta,
        ingff,
        outgb)
    return cmd


if __name__ == "__main__":  # pragma: no cover
    args = get_args()
    # allow user to give relative paths
    output_root = os.path.abspath(os.path.expanduser(args.output))
    try:
        os.makedirs(output_root, exist_ok=False)
    except OSError:
        print("Output directory already exists; exiting...")
        sys.exit(1)
    t0 = time.time()
    log_path = os.path.join(output_root,
                            str("{0}_riboScan_log.txt".format(
                                time.strftime("%Y%m%d%H%M"))))
    logger = set_up_logging(verbosity=args.verbosity,
                            outfile=log_path,
                            name=__name__)
    # # log version of riboSeed, commandline options, and all settings

    logger.info("Usage:\n{0}\n".format(" ".join([x for x in sys.argv])))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))

    output_file = os.path.join(
        output_root,
        "scannedScaffolds.gb")
    ##  get and check list of input files
    fastas = glob.glob(str(args.contigs_dir + "*" + args.ext))
    if len(fastas) == 0:
        logger.error("No fasta files in %s with extention %s! Exiting",
                     args.contigs_dir, args.ext)
        sys.exit(1)

    gb_list = []
    for fasta in fastas:
        with open(fasta, 'r') as f:
            header = f.readline().strip()
        accession = parse_fasta_header(header)
        logger.info("Accession for %s: %s", fasta, accession)
        barrnap_cmd = make_barrnap_cmd(
            infasta=fasta,
            outgff=os.path.join(output_root, "{0}.gff".format(accession)),
            exe=args.barrnap_exe,
            thresh=args.id_thresh,
            kingdom=args.kingdom)
        logger.info("running barrnap cmd: %s", barrnap_cmd)
        subprocess.run(barrnap_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
        logger.info("grooming gff")
        tagged_gff = add_locus_tags_to_gff(
            gff=os.path.join(output_root, "{0}.gff".format(accession)),
            acc=accession)
        unfinished_gb = os.path.join(output_root, "{0}_pre.gb".format(accession))
        seqret_cmd = make_seqret_cmd(exe=args.seqret_exe,
                                     outgb=os.path.join(output_root, "{0}_pre.gb".format(accession)),
                                     ingff=tagged_gff, infasta=fasta)
        subprocess.run(seqret_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
        append_accession_and_version(accession=accession,
                                     ingb=unfinished_gb,
                                     finalgb=os.path.join(output_root, "{0}.gb".format(accession)))
        gb_list.append(os.path.join(output_root, "{0}.gb".format(accession)))
    ### no more fastas
    combine_gbs(finalgb=output_file, gb_list=gb_list)
    # Report that we've finished
    logger.info("Done: %s", time.asctime())
    logger.info("combined scaffolds can be found here: %s", output_file)
    logger.info("Time taken: %.3fs" % (time.time() - t0))
