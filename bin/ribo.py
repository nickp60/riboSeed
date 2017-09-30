#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import sys
# print("\n".join(sys.path))
from riboSeed import __version__
from riboSeed.scripts import run_riboSeed
from riboSeed import riboScan
from riboSeed import riboSelect
from riboSeed import riboSeed as seed
from riboSeed import riboScore
from riboSeed import riboSketch
from riboSeed import riboSim

helpmsg = [
    "riboSeed v" + __version__,
    "Contact: Nick Waters <nickp60@gmail.com>",
    "",
    "Usage:   ribo <command> [options]",
    "",
    "Availible commands:",
    "-run: execute the whole pipeline the whole pipeline",
    "-scan: reannotate rDNAs in a FASTA file isong help ",
    "-select: group rDNA annotations into operons",
    "etc"
]


def main():
    if len(sys.argv) == 1:
        print("\n".join(helpmsg))
        sys.exit(1)
    if "-v" in sys.argv[1]:
        print("version")
        sys.exit(1)
    modules_dict = {
        "run": run_riboSeed,
        "scan": riboScan,
        "select": riboSelect,
        "seed": seed,
        "sketch": riboSketch,
        "score": riboScore,
        "sim": riboSim}
    if sys.argv[1] not in modules_dict.keys():
        print("help")
        sys.exit(1)
    thisprogram = sys.argv[1]
    print(thisprogram)
    args = modules_dict[thisprogram].get_args()
    sys.exit(modules_dict[thisprogram].main(args))


if __name__ == '__main__':
    main()
