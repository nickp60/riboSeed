#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import sys
# print("\n".join(sys.path))
import importlib
import riboSeed

helpmsg = [
    "riboSeed v" + riboSeed.__version__,
    "Contact: Nick Waters <nickp60@gmail.com>",
    "",
    "Usage:   ribo <command> [options]",
    "",
    "Availible commands:",
    "-run: execute the whole pipeline the whole pipeline",
    "-scan: reannotate rRNAs in a FASTA file isong help ",
    "-select: group rRNA annotations into rDNA operons",
    "-seed: perform de fere novo assembly",
    "-snag: extract rDNA regions and plot entropy",
    "-sim: perform simulations used in manuscript",
    "-sketch: plot results from a de fere novo assembly",
    "-score: score batches of assemblies with BLASTn",
    "-swap: swap contigs from assemblies"
]


def main():
    assert ((sys.version_info[0] == 3) and (sys.version_info[1] >= 5)), \
        "Must use python3.5 or higher!"
    if len(sys.argv) == 1:
        print("\n".join(helpmsg))
        sys.exit(1)
    if "-v" in sys.argv[1]:
        print("version")
        sys.exit(1)
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
    if sys.argv[1] not in modules_dict.keys():
        print(sys.argv[1] + " is not an available program. See the list below:")
        print("\n".join(helpmsg))
        sys.exit(1)
    modname = sys.argv[1]
    this_module = importlib.import_module("riboSeed." + modules_dict[modname])

    args = this_module.get_args()
    sys.exit(this_module.main(args))


if __name__ == '__main__':
    main()
