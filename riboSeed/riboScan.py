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

# need this line for unittesting
sys.path.append(os.path.join('..', 'riboSeed'))
from pyutilsnrw.utils3_5 import set_up_logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# --------------------------- classes --------------------------- #


# --------------------------- methods --------------------------- #


def get_args():  # pragma: no cover
    """#TODO:     for cli mods:
    http://stackoverflow.com/questions/18025646/
         python-argparse-conditional-requirements
    make this able to handle different library types such as two unpaired runs
    """
    parser = argparse.ArgumentParser(
        description="Given a directory of one or more chromosomes as fasta " +
        "files, this facilitates reannotation of rDNA regions with Barrnap " +
        " and outputs all sequences as a single, annotated genbank file",
        add_help=False)  # to allow for custom help
    parser.add_argument("contigs", action="store",
                        help="either a (multi)fasta or a directory " +
                        "containing one or more chromosomal " +
                        "sequences in fasta format")

    # taking a hint from http://stackoverflow.com/questions/24180527
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)

    # had to make this faux "optional" parse so that the named required ones
    # above get listed first
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-e", "--extension", dest='ext',
                          action="store",
                          help="extension of the chromosomal sequences, " +
                          "usually '.fasta' or similar",
                          default=".fa",
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
    optional.add_argument("-b", "--barrnap_exe", dest='barrnap_exe',
                          action="store",
                          help="path to barrnap executable; " +
                          "default: %(default)s", default="barrnap",
                          type=str)
    optional.add_argument("-n", "--name", dest='name',
                          action="store",
                          # default="contig",
                          help="name to give the contig files; "
                          "default: infered from file",
                          type=str)
    optional.add_argument("-c", "--cores", dest='cores',
                          action="store",
                          help="number of threads/cores to use; " +
                          "default: %(default)s", default=2,
                          choices=[1, 2, 4, 8, 16], type=int)
    optional.add_argument("-s", "--seqret_exe", dest='seqret_exe',
                          action="store",
                          help="path to seqret executable, usually " +
                          "installed with emboss; " +
                          "default: %(default)s", default="seqret",
                          type=str)
    optional.add_argument("-m", "--min_length", dest='min_length',
                          action="store",
                          help="skip the scaffold if its shorter than this " +
                          "default: %(default)s", default=0,
                          type=int)
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
        try:  # this handles normal flavored fastas with '>seqid descr'
            accession = re.search(r">(.*?) .*", first_line).groups()[0]
        except AttributeError:  # this handles headers with only '>seqid\n'
            accession = first_line.strip().replace(">", "")
    if accession is None:
        raise ValueError(
            "unable to extract accession from first line in file: \n%s".format(
                first_line))
    return accession


def make_barrnap_cmd(infasta, outgff, exe, thresh, kingdom, threads=1):
    assert shutil.which(exe) is not None, "barrnap executable not found!"
    assert thresh > 0 and thresh < 1, "Thresh must be between 0 and 1!"
    if exe.endswith("py"):
        # ensure running python barrnap uses >3.5
        pyexe = str(sys.executable + " ")
    else:
        pyexe = ""
    cmd = "{0}{1} -k {2} {3} --reject {4} --threads {5} > {6}".format(
        pyexe,
        shutil.which(exe),
        kingdom,
        infasta,
        thresh,
        threads,
        outgff
        )
    return cmd


def add_locus_tags_to_gff(gff, acc):
    LOCUS = 0
    gff_list = []
    with open(gff, 'r') as g:
        for idx, line in enumerate(g):
            line = line.strip()
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
            # print(str(l))
            outgff.write(str(l + "\n"))
    return(new_gff)


def combine_gbs(finalgb, gb_list):
    with open(finalgb, 'w') as outfile:
        for idx, fname in enumerate(gb_list):
            with open(fname) as infile:
                for line in infile:
                    if line.startswith("##") and idx != 0:
                       continue
                    outfile.write(line)


def append_accession_and_version(accession, ingb, finalgb):
    """updated to add a dummy GI number as well
    """
    with open(finalgb, 'w') as outfile:
        with open(ingb) as infile:
            for idx, line in enumerate(infile):
                if idx == 1:
                    repline = "ACCESSION   {0}\nVERSION     {0}{1}{2}".format(
                        accession,
                        "  GI:0000000\n",
                        line)
                    outfile.write(repline)
                else:
                    outfile.write(line)


def make_seqret_cmd(exe, outgb, ingff, infasta):
    assert shutil.which(exe) is not None, "seqret executable not found!"
    cmd = str(
        "{0} -sequence {1} -feature -fformat gff3 -fopenfile {2} " +
        "-osformat genbank -auto -outseq {3}"
    ).format(
        shutil.which(exe),
        infasta,
        ingff,
        outgb)
    return cmd


def checkSingleFasta(fasta, logger):
    with open(fasta, 'r') as f:
        counter = 0
        for rec in SeqIO.parse(f, "fasta"):
            counter = counter + 1
            if counter > 1:
                logger.error(
                    "Fasta has multiple entries! " +
                    "use splitMultifasta.sh to create single-entry " +
                    "fastas from it. exiting...")
                sys.exit(1)


def getFastas(inp, output_root, ext, name, logger):
    if not os.path.isdir(os.path.expanduser(inp)):
        if not os.path.isfile(os.path.expanduser(inp)):
            logger.error("'%s' is not a valid directory or file!",
                         os.path.expanduser(inp))
            sys.exit(1)
        else:
            if not os.path.getsize(os.path.expanduser(inp)) > 0:
                logger.error("'%s' looks like an empty file!",
                             os.path.expanduser(inp))
                sys.exit(1)
            logger.info("spliting multifasta into multiple fastas " +
                        "for easier processing")
            splitMultifasta(multi=inp, output=output_root, name=name,
                            logger=logger)

            fastas = glob.glob(os.path.join(output_root, "contigs",
                                            "*" + ext))
    else:
        fastas = glob.glob(os.path.join(os.path.expanduser(inp),
                                        "*" + ext))
    if len(fastas) == 0:
        logger.error("No fasta files in %s with extention %s! Exiting",
                     inp, ext)
        sys.exit(1)
    # unify the output of multifasta
    combine_gbs(finalgb=os.path.join(output_root, "scannedScaffolds.fa"),
                gb_list=fastas)
    return(fastas)


def splitMultifasta(multi, output, name, dirname="contigs", logger=None):
    """regex stolen from SO
    name is the name of the file, output is the parent dir for the output dir
    """
    idlist = []
    assert logger is not None, "must use logging"
    os.makedirs(os.path.join(output, "contigs"))
    with open(os.path.expanduser(multi), "r") as mf:
        for idx, rec in enumerate(SeqIO.parse(mf, "fasta")):
            logger.debug("record %d: %s", idx, rec.id)
            # logger.debug(rec)
            if name is not None:
                fname = name
            else:
                fname = rec.id
            # if fname in idlist:
            fname = fname + "_" + str(idx)
            if not re.match("^[a-zA-Z0-9_\.]*$", fname):
                logger.error("Problem with header %s", fname)
                logger.error("file header contains special characters! " +
                             "Valid characters are in set [a-zA-Z0-9_]; " +
                             "rename with the --name option to prevent " +
                             "later issues. I'm sorry, its a pain...\n" +
                             "Exiting...")
                sys.exit(1)
            with open(os.path.join(output, "contigs",
                                   fname + ".fa"), "w") as outf:
                renamed_rec = SeqRecord(rec.seq, id=fname, description="")
                SeqIO.write(renamed_rec, outf, "fasta")
            idlist.append(fname)
    logger.debug("rewrote the following records")
    logger.debug("\n".join(idlist))

if __name__ == "__main__":
    args = get_args()
    # allow user to give relative paths
    output_root = os.path.abspath(os.path.expanduser(args.output))
    try:
        os.makedirs(output_root, exist_ok=False)
    except OSError:
        print("Output directory already exists; exiting...")
        sys.exit(1)
    t0 = time.time()
    log_path = os.path.join(output_root, "riboScan.log")
    logger = set_up_logging(verbosity=args.verbosity,
                            outfile=log_path,
                            name=__name__)
    # # log version of riboSeed, commandline options, and all settings

    logger.info("Usage:\n{0}\n".format(" ".join([x for x in sys.argv])))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))

    output_file = os.path.join(output_root, "scannedScaffolds.gb")
    output_file_gff = os.path.join(output_root, "scannedScaffolds.gff")
    ##  get and check list of input files
    fastas = getFastas(inp=args.contigs, output_root=output_root,
                       ext=args.ext, name=args.name, logger=logger)
    # if len(fastas) == 0:
    #     logger.error("No fasta files in %s with extention %s! Exiting",
    #                  args.contigs, args.ext)
    #     sys.exit(1)

    gb_list = []
    gff_list = []  # for tagged gffs, that is
    for fasta in sorted(fastas):
        # check for multiple entry fastas.
        # This can still happen if mutlifastas are in dir input
        checkSingleFasta(fasta, logger=logger)
        with open(fasta, 'r') as f:
            rec = next(SeqIO.parse(f, "fasta"))
        accession = rec.id
        if args.min_length > len(rec.seq):
            logger.debug("Skipping %s: must be %d bp, set by --min_length",
                         accession, args.min_length)
            continue
        logger.info("Accession for %s: %s", fasta, accession)
        barrnap_cmd = make_barrnap_cmd(
            infasta=fasta,
            outgff=os.path.join(output_root, "{0}.gff".format(accession)),
            exe=args.barrnap_exe,
            threads=args.cores,
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
        unfinished_gb = os.path.join(
            output_root, "{0}_pre.gb".format(accession))
        seqret_cmd = make_seqret_cmd(
            exe=args.seqret_exe,
            outgb=os.path.join(output_root, "{0}_pre.gb".format(accession)),
            ingff=tagged_gff, infasta=fasta)
        subprocess.run(seqret_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
        append_accession_and_version(
            accession=accession,
            ingb=unfinished_gb,
            finalgb=os.path.join(output_root, "{0}.gb".format(accession)))
        gb_list.append(os.path.join(output_root, "{0}.gb".format(accession)))
        gff_list.append(tagged_gff)
    ### no more fastas
    combine_gbs(finalgb=output_file, gb_list=gb_list)
    # this is why you give things intellegent method names
    combine_gbs(finalgb=output_file_gff, gb_list=gff_list)
    # Report that we've finished
    logger.info("Done: %s", time.asctime())
    logger.info("combined scaffolds can be found here: %s", output_file)
    logger.info("Time taken: %.3fs" % (time.time() - t0))
