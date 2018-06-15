#!/bin/bash
set -euo


if [ -z "$3" ]
then
    echo "3 required arguments: forward reads, reverse reads, name. \n USAGE: mini_assembly.sh reads1 reads2 NAME"
    exit 1
else
    NAME="$3"
fi



OUTDIRBASE="./"
OUTDIR="${OUTDIRBASE}`date +%F`_mini_assembly_${NAME}/"
mkdir ${OUTDIR}
mkdir ${OUTDIR}/tmp/

average_nucleotide_identity.py -i ./subgenomes/ -g -o pyani_test -v

col=$(python ~/GitHub/riboSeed/scripts/colnum.py ./ani_results.tab )
target=$(head -n 1 ./pyani_test/ANIm_percentage_identity.tab | cut -f "$(($col+1))")
