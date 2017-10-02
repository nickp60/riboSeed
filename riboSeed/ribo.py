#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import sys
# print("\n".join(sys.path))
import importlib
import riboSeed

CITATION = ""
helpmsg = [
    "riboSeed v" + riboSeed.__version__,
    "Contact: Nick Waters <nickp60@gmail.com>",
    "Description: A suite of tools to perform de fere novo assembly to bridge",
    "             gaps caused by rDNA repeats",
    "",
    "Usage:   ribo <command> [options]",
    "",
    "Available commands:",
    " -run      execute the whole pipeline (scan, select, seed, sketch, and score)",
    " -scan     reannotate rRNAs in a FASTA file isong help ",
    " -select   group rRNA annotations into rDNA operons",
    " -seed     perform de fere novo assembly",
    " -snag     extract rDNA regions and plot entropy",
    " -sim      perform simulations used in manuscript",
    " -sketch   plot results from a de fere novo assembly",
    " -score    score batches of assemblies with BLASTn",
    " -swap     swap contigs from assemblies"
]


def main(args):
    assert ((sys.version_info[0] == 3) and (sys.version_info[1] >= 5)), \
        "Must use python3.5 or higher!"
    if len(args) == 1:
        print("\n".join(helpmsg))
        return 0
    if "-v" in args[1]:
        print(riboSeed.__version__)
        return 0
    modules_dict = {
        "run": "run_riboSeed",
        "scan": "riboScan",
        "select": "riboSelect",
        "snag": "riboSnag",
        "seed": "riboSeed",
        "sketch": "riboSketch",
        "score": "riboScore",
        "sim": "riboSim",
        "swap": "riboSwap",
        "config": "make_riboSeed_config"}
    if args[1] not in modules_dict.keys():
        print("Error:" + args[1] +
              " is not an available program. See the list below:\n")
        print("\n".join(helpmsg))
        return 1
    modname = args[1]
    this_module = importlib.import_module("riboSeed." + modules_dict[modname])

    these_args = this_module.get_args()
    return this_module.main(these_args)


if __name__ == '__main__':
    args = sys.argv
    main(args)