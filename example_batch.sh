#!/bin/bash
# version 0.0.1

# Runs a riboseed assembly with default settings

# output scanScaffolds_combined.gb in current directory
echo 'USAGE: /path/to/genome.gb path/to/genome.fasta path/to/read1 path/to/read2 /path/to/outdir/ iterations'
echo 'example: $ barrnap'
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3"] || [ -z "$4"]|| [ -z "$5"] ||[ -z "$6"] ||[ -z "$7"]
then
    echo "mandatory arguments: genome.fasta,genome.gb, read1, read2, output_dir, iterations, and flanking_width"
    exit 1
fi

if [ -d "$5" ]; then
    echo "output dir exists!"
    echo "play things safe, make a new directory for results"
    exit 1
fi
GB="$1"
FASTA="$2"
READ1="$3"
READ2="$4"
OUTDIR="$5"
ITERATIONS="$6"
FLANK="$7"
## make  dirs
mkdir ${OUTDIR}

## RiboSelect
riboSelect.py ${GB} -o ${OUTDIR}select/
## RiboSnag
riboSnag.py ${GB} ${OUTDIR}select/
# RiboSeed
riboSeed.py ./pseudata/snag/ -F ./toy_reads1.fq -R toy_reads2.fq -r ./data/littlegenome.fa -v 1 -i 3 --ref_as_contig trusted -c 4 -o ./pseudata/seed7777/
