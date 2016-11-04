#!/bin/bash
# version 0.0.4
#$ -cwd
#$ -j yes
#$ -V

echo "Hostname: $HOSTNAME"
echo "Task started `date`"

set -e # exit on error!
# Runs a riboseed assembly with default settings

# output scanScaffolds_combined.gb in current directory
echo 'USAGE: /path/to/genome.gb path/to/genome.fasta path/to/read1 path/to/read2 /path/to/outdir/ n_iterations n_flanking n_cores'

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] || [ -z "$8" ] ;
then
    echo "All mandatory arguments: genome.gb, genome.fasta, read1, read2, output_dir, iterations, flanking_width, virtenve_exe"
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
VENVEXE="$8"
## make  dirs
mkdir ${OUTDIR}
## enter virtual environment
source ${VENVEXE}
## RiboSelect
python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ${GB} -o ${OUTDIR}select/
## RiboSnag
python3.5 ~/GitHub/riboSeed/riboSeed/riboSnag.py ${GB} ${OUTDIR}select/riboSelect_grouped_loci.txt -o ${OUTDIR}snag/ -l $FLANK
# RiboSeed
python3.5 ~/GitHub/riboSeed/riboSeed/riboSeed.py ${OUTDIR}snag/ -F ${READ1} -R ${READ2} -r ${FASTA} -v 1 -i ${ITERATIONS} -o ${OUTDIR}seed/
