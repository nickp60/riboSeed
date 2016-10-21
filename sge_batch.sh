#!/bin/bash

#$ -cwd
#$ -j yes
#$ -V

echo "Hostname: $HOSTNAME"
echo "Task started `date`"

GB="$1"
FASTA="$2"
READ1="$3"
READ2="$4"
OUTDIR="$5"
ITERATIONS="$6"
FLANK="$7"
CORES="$8"

# include any commands and/or other
# bash scripting instructions here
./example_batch.sh $GB $FASTA $READ1 $READ2 $OUTDIR $ITERATIONS $FLANK $CORES

echo "Task ended `date`"
