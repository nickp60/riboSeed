#!/bin/bash
# version 0.0.1
# plagerized awk bit from https://gist.github.com/astatham/621901
# argument 1 is path to a multifasta file
# takes a multifasta, splits into multifasta

echo 'USAGE: /path/to/contigs.fasta'

# replace extraneous stuff with sed, run awk
cat "$1" | sed 's/ <unknown description>//' | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}
        print $0 > filename
}'
