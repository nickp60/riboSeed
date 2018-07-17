#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# Copyright 2017, National University of Ireland and The James Hutton Insitute
# Author: Nicholas Waters
#
# This code is part of the riboSeed package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

import argparse
import sys
import time
import os
import shutil
import pkg_resources


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="Generate a configuration yaml file for your riboSeed run",
        add_help=True)
    parser.add_argument("-o", "--outdir", dest='outdir', action="store",
                        help="output directory; " +
                        "default: %(default)s", default=os.getcwd(),
                        type=str, required=False)
    parser.add_argument("-n", "--name", dest='name', action="store",
                        help="name of config file; " +
                        "default: timestamped", default=None,
                        type=str, required=False)
    args = parser.parse_args(sys.argv[2:])
    return args


def make_config_header():
    hashbang = []
    copy_header = [
        "# Copyright 2017, National University of Ireland and The James Hutton Insitute",
        "# Author: Nicholas Waters",
        "#",
        "# This code is part of the riboSeed package, and is governed by its licence.",
        "# Please see the LICENSE file that should have been included as part of",
        "# this package.\n\n"]
    t0 = time.asctime()
    noteline = "# NOTE: 'null' will be interpretted as a python None "
    timestamp = "# This config file was generated " + str(t0) + "\n\n"
    hashbang.extend(copy_header)
    hashbang.append(noteline)
    hashbang.append(timestamp)
    return(hashbang)


def config_exes():
    # find all needed system requiremnts
    req_programs = [("BARRNAP_EXE", "barrnap"),
                    ("SEQRET_EXE", "seqret"),
                    ("SPADES_EXE", "spades.py"),
                    ("BWA_EXE", "bwa"),
                    ("SAMTOOLS_EXE", "samtools"),
                    ("BLAST_EXE", "blastn"),
                    ("MAKEBLASTDB_EXE", "makeblastdb")]
    config_lines = [
        "#------------------------#",
        "##  Required programs   ##",
        "#------------------------#"]
    for k, v in req_programs:
        config_lines.append("# executable for " + v)
        if shutil.which(v):
            config_lines.append(k + " : " + shutil.which(v) + "\n")
        else:
            sys.stderr.write("Warning! No executable found for " + v + "; " +
                             "Please add it to the config file manually\n")
            config_lines.append(k + " : null\n" )

    # # find all optional sys requirements
    opt_programs = [("QUAST_EXE", "quast"),
                    ("SMALT_EXE", "smalt"),
                    ("PRANK_EXE", "prank"),
                    ("MAFFT_EXE", "mafft"),
                    ("BCFTOOLS_EXE", "bcftools")]
    config_lines.extend([
        "#------------------------#",
        "##  Optional programs   ##",
        "#------------------------#"])
    for k, v in opt_programs:
        config_lines.append("# executable for " + v)
        if shutil.which(v):
            config_lines.append(k + " : " + shutil.which(v) + "\n")
        else:
            config_lines.append(k + " : null\n" )
    config_lines.extend(config_mauve())
    return config_lines


def config_mauve():
    mauve_config_lines = [
        "#------------------------#",
        "##   Mauve (Optional)   ##",
        "#------------------------#",
        "# Muave is used by riboSketch for orienting contigs and ",
        "# determining synteny.  riboSketch uses the Mauve.jar java program,",
        "# but usually only mauveAligner (in the platform specific subdir)",
        "# is actually added to the path after installation."
    ]

    mauve_aligner_name = "mauveAligner"
    if shutil.which(mauve_aligner_name):
        mauve_aligner_exe = shutil.which(mauve_aligner_name)
        # This assumes Mauve keeps installing the exe' in a platform specific
        # subfolder, and leaves the jar in the parent of that platform folder
        mauve_jar =  os.path.join(
            os.path.dirname(os.path.dirname(shutil.which("mauveAligner"))),
            "Mauve.jar")
        if not os.path.exists(mauve_jar):
            mauve_jar=None
    else:
        mauve_aligner_exe = None
        mauve_jar = None

    programs = [("MAUVE_ALIGNER_EXE", mauve_aligner_exe, "mauveAligner"),
                ("MAUVE_JAR", mauve_jar, "Mauve.jar")]

    for k, v, exe in programs:
        mauve_config_lines.append("# executable for " + exe)
        if v is not None:
            mauve_config_lines.append(k + " : " + v + "\n")
        else:
            mauve_config_lines.append(k + " : null\n" )
    return mauve_config_lines


def config_scan_defaults():
    scan_params = [
        ("SCAN_ID_THRESH", "0.5",
         "--id_thresh: partial rRNA hits below this threshold will be ignored"),
        ("SCAN_CONTIG_NAME", "null",
         "--name: stem to give the files generated by riboScan  "),
        ("SCAN_MIN_LENGTH", "0",
         "--min_length: skip annotating if contig shorter than this"),
        ("SCAN_VERBOSITY", "2",
         "-v: verbosity for riboScan")]
    scan_lines = [
        "#------------------------#",
        "##  riboScan Parameters ##",
        "#------------------------#"]
    for k, v, h in scan_params:
        scan_lines.append("# " + h)
        scan_lines.append(k + " : " + v + "\n")
    return scan_lines


def config_select_defaults():
    select_params = [
        ("SELECT_FEATURE", "'rRNA'",
         "--feature: which annotations to pay attention to; \n# " +
         "barrnap uses 'rRNA', but others may use 'RRNA', etc"),
        ("SELECT_VERBOSITY", "2",
         "-v: verbosity for riboSelect")]
    select_lines = [
        "#------------------------#",
        "## riboSelect Parameters #",
        "#------------------------#"]
    for k, v, h in select_params:
        select_lines.append("# " + h)
        select_lines.append(k + " : " + v + "\n")
    return select_lines


def config_snag_defaults():
    snag_params = [
        ("SNAG_NAME", "'name'",
         "name to use when plotting the results"),
        ("SNAG_FLANKING", "'name'",
         "length (bps) of flanking regions to analyze"),
        ("SNAG_MSA_KMERS", "false",
         ""),
        ("SNAG_SKIP_KMERS", "false",
         ""),
        ("SNAG_SKIP_BLAST", "false",
         ""),
        ("SNAG_LINEAR", "false",
         ""),
        ("SNAG_TITLE", "title",
         ""),
        ("SNAG_NO_REVCOMP", "false",
         ""),
        ("SNAG_JUST_EXTRACT", "false",
         ""),
        ("SNAG_MSA_TOOL", "mafft",
         ""),
        ("SNAG_PADDING", "5000",
         ""),
        ("SNAG_VERBOSITY", "2",
         "-v: verbosity for riboSnag")]
    snag_lines = [
        "#------------------------#",
        "## riboSnag Parameters #",
        "#------------------------#",
        "# TODO:the documentation here could be improved"]
    for k, v, h in snag_params:
        snag_lines.append("# " + h)
        snag_lines.append(k + " : " + v + "\n")
    return snag_lines

def config_seed_defaults():
    seed_params = [
        ("SEED_MAP_METHOD", "bwa",
         "--method_for_map: which mapper to use, (default bwa)"),
        # ("SEED_SCORE_MIN", "null",
        #  "--score_min: If using smalt, this sets the '-m' param; \n# " +
        #  "default with smalt is inferred from \n# " +
        #  "read length. If using BWA, reads mapping with AS\n# " +
        #  "score lower than this will be rejected\n# " +
        #  "; default with BWA is half of read length"),
        ("SEED_MIN_ASSEMBLY_LENGTH", "6000",
         "if initial SPAdes assembly largest contig \n# " +
         "is not at least as long as --min_assembly_len, \n# " +
         "reject. Set this to the length of the seed \n# " +
         "sequence; if it is not achieved, seeding across \n# " +
         "regions will likely fail; default: %(default)s"),
        ("SEED_INCLUDE_SHORTS", "false",
         "if assembled contig is smaller than --min_assembly_len, contig\n#" +
         "will still be included in assembly; default: inferred"),
        ("SEED_SUBTRACT", "false",
         "if --subtract reads already used in previous\n# " +
         "round of subassembly will not be included in \n# " +
         "subsequent rounds.  This can lead to problems \n# " +
         "with multiple mapping and inflated coverage."),
        ("SEED_SKIP_CONTROL", "false",
         "if --skip_control, no de novo \n# " +
         "assembly will be done; default: %(default)s"),
        # ("SEED_REF_AS_CONTIG", "infer",
        #  "ignore: reference will not be used in \n# " +
        #  "subassembly. trusted: SPAdes will use the seed\n# " +
        #  " sequences as a --trusted-contig; untrusted: \n# " +
        #  "SPAdes will treat as --untrusted-contig. \n# " +
        #  "infer: if mapping \n#" +
        #  "percentage over 80%: 'trusted', else 'untrusted'\n# " +
        #  "See SPAdes docs for details.  default: infer"),
        ("SEED_TARGET_LEN", "null",
         "If set, iterations will continue until \n# " +
         "contigs reach this length, or max iterations (\n#" +
         "set by --iterations) have been completed. Set as \n#" +
         "fraction of original seed length by giving a \n#" +
         "decimal between 0 and 5, or set as an absolute \n#" +
         "number of base pairs by giving an integer greater\n#" +
         " than 50. Not used by default"),
        ("SEED_MAPPER_ARGS", "-L 0,0 -U 0 -a",
         "Submit custom parameters to mapper. \n#" +
         "And by mapper, I mean bwa, cause we dont support \n#" +
         "this option for SMALT, sorry. \n#" +
         "This requires knowledge of your chosen mapper's \n#" +
         "optional arguments. Proceed with caution!"),
        ("SEED_SMALT_SCORING", "match=1,subst=-4,gapopen=-4,gapext=-3",
         "If mapping with SMALT, \n#" +
         "submit custom smalt scoring via smalt -S \n#" +
         "scorespec option"),
        ("SEED_VERBOSITY", "2",
         "-v: verbosity for riboSeed"),
        ("SEED_INITIAL_CONSENSUS", "False",
         "--initial_consensus enables a consensus approach to the first\n#" +
         " iterations subassemblies, rather than a de Druijn graph de non\n#" +
         " assembly with SPAdes")]
    seed_lines = [
        "#------------------------#",
        "##  riboSeed Parameters ##",
        "#------------------------#"]
    for k, v, h in seed_params:
        seed_lines.append("# \n# " + h)
        seed_lines.append(k + " : " + v + "\n")
    return seed_lines


def config_sketch_defaults():
    sketch_params = [
        ("SKETCH_ASSEMBLY_EXT", "'fasta'",
         "-f extension of assemblies, usually fasta"),
        ("SKETCH_REF_EXT", "'gb'",
         "-g extension of reference, usually gb"),
        ("SKETCH_VERBOSITY", "2",
         "-v: verbosity for riboSketch")]
    sketch_lines = [
        "#------------------------#",
        "## riboSketch Parameters #",
        "#------------------------#"]
    for k, v, h in sketch_params:
        sketch_lines.append("# \n# " + h)
        sketch_lines.append(k + " : " + v + "\n")
    return sketch_lines


def config_score_defaults():
    score_lines = [
        "#------------------------#",
        "## riboScore Parameters ##",
        "#------------------------#",
        "# NOTE: other than verbosity, there are no parameters that are ",
        "# unique to riboScore, as all the arguments have den defined in ",
        "# previous programs above. See run_riboSeed.py for more info as to ",
        "# which args are used."]
    score_lines.append("# -v: verbosity for riboSketch")
    score_lines.append("SCORE_VERBOSITY : 2\n")
    return score_lines

def config_stack_defaults():
    stack_params = [
        ("STACK_N_SAMPLES", "10",
         "number of samples to compare rRNA region depth to"),
        ("STACK_INFER", "false",
         "how to handle naming of regions"),
        ("STACK_VERBOSITY", "2",
         "-v: verbosity for riboSketch")]
    stack_lines = [
        "#------------------------#",
        "## riboStack Parameters #",
        "#------------------------#"]
    for k, v, h in stack_params:
        stack_lines.append("# \n# " + h)
        stack_lines.append(k + " : " + v + "\n")
    return stack_lines

def config_spec_defaults():
    spec_params = [
        ("SPEC_MIN_CONTIG_LEN", "75",
         "minimum length of contig to consider for depth calculations"),
        ("SPEC_MIN_ANCHOR_LENGTH", "500",
         "The minimum length for an 'anchor contig'; see riboSpec docs"),
        ("SPEC_THRESHOLD", "1000",
         "Threshold for minimun length of valid graph paths away from rDNA "),
        ("SPEC_MEDIUM_LENGTH_THRESHOLD", "400",
         "Threshold for defining medium-length contigs"),
        ("SPEC_BARRNAP_LENGTH_THRESHOLD", ".75",
         "Threshold barrnap uses for length of hits"),
        ("SPEC_VERBOSITY", "2",
         "-v: verbosity for riboSketch")]
    spec_lines = [
        "#------------------------#",
        "## riboSpec Parameters #",
        "#------------------------#"]
    for k, v, h in spec_params:
        spec_lines.append("# \n# " + h)
        spec_lines.append(k + " : " + v + "\n")
    return spec_lines


def config_run_defaults():
    run_lines = [
        "#------------------------#",
        "##    Run Parameters    ##",
        "#------------------------#",
        "# NOTE: This section will only be filled in when run by run_riboSeed",
        "# as this contains info specific to a particular run.  This feature ",
        "# is aimed at improving reproducibility.  See the docs for details ",
        "# about the args"]
    return run_lines


def write_config(header, outfile):
    with open(outfile, "w") as of:
        for hline in header:
            of.write(hline)
            of.write("\n")


def append_config(lines, outfile):
    with open(outfile, "a") as of:
        for line in lines:
            of.write(line)
            of.write("\n")


def main(args):
    header = make_config_header()
    if args.name is None:
        outfile = os.path.join(
            os.path.abspath(os.path.expanduser(args.outdir)),
            str(time.strftime("%Y-%m-%dT%H%M") +
                "_riboSeed_config.yaml"))
    else:
        outfile = os.path.join(
            os.path.abspath(os.path.expanduser(args.outdir)),
            str(args.name + ".yaml"))
    write_config(header=header, outfile=outfile)
    #exes
    append_config(lines=config_exes(), outfile=outfile)
    # riboScan
    append_config(lines=config_scan_defaults(), outfile=outfile)
    # riboSelect
    append_config(lines=config_select_defaults(), outfile=outfile)
    # riboSeed
    append_config(lines=config_seed_defaults(), outfile=outfile)
    # riboSketch
    append_config(lines=config_sketch_defaults(), outfile=outfile)
    # riboScore
    append_config(lines=config_score_defaults(), outfile=outfile)
    # riboStack
    append_config(lines=config_stack_defaults(), outfile=outfile)
    # riboSpec
    append_config(lines=config_spec_defaults(), outfile=outfile)
    # riboSnag
    append_config(lines=config_snag_defaults(), outfile=outfile)
    # run_riboSeed
    append_config(lines=config_run_defaults(), outfile=outfile)
    # this lines allows run_riboseed to find the path to the new config
    sys.stdout.write(outfile)
    # this line allows us to just grab the main function in run_riboSeed
    return outfile
