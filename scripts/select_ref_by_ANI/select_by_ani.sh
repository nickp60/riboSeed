#!/bin/bash
set -euo


if [ -z "$1" ]
then
    echo "1 required argument: \n USAGE: selecy_by_ani.sh genomes_dir"
    exit 1
else
    GENOMESDIR="$1"
fi



OUTDIRBASE="./"
OUTDIR="${OUTDIRBASE}`date +%F`_ANI/"
mkdir ${OUTDIR}

average_nucleotide_identity.py -i $GENOMESDIR -g -o $OUTDIR/pyani

python ~/GitHub/riboSeed/scripts/select_ref_by_ANI/colnum.py $OUTDIR/pyani/ANIm_percentage_identity.tab

# target=$(head -n 1 $OUTDIR/tmp/pyani/ANIm_percentage_identity.tab | cut -f "$(($col+1))")
# target=$(tail -n +2 $OUTDIR/tmp/pyani/ANIm_percentage_identity.tab | awk '{print $(($col+1)) "\t" $1}' | sort -n -r | tail -n +2 | head -n 1 | cut -f 2
# echo $target
