#!/bin/bash
# version 0.0.2
# plagerized awk bit from https://gist.github.com/astatham/621901
# argument 1 is path to a multifasta file
# takes a multifasta, splits into multiple fasta files

echo 'USAGE: /path/to/contigs.fasta'

# replace extraneous stuff with sed, run awk
cat "$1" | awk '{print $1}' | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0, 2) ".fa")}
        print $0 > filename
}'
