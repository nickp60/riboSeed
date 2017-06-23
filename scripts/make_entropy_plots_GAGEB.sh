#!/bin/bash
# version 0.2.0
IFS=$'\n\t'

# process original
# only need to do this once;
# we arent doing indels so the coords are the same
OUTDIRBASE="./"
OUTDIR="${OUTDIRBASE}`date +%F`_GAGEB_entropy"

mkdir ${OUTDIR}
mkdir "${OUTDIR}/figures/"
accessions=`cat ./scripts/genomes_accessions.txt`
filecounter=10
for line in ${accessions}
do
    IFS=$'\n\t '
    mkdir ${OUTDIR}/${filecounter}_dir
    for i in $line
    do
	# for each accession in line, if more than one
	if [ -e "${OUTDIR}/${filecounter}/${i}.fasta" ]
	then
	    echo "using local copy of reference"
	else
	    get_genomes.py -q $i -o ${OUTDIR}/${filecounter}_dir
	fi
    done

    python3.5 ~/GitHub/riboSeed/riboSeed/riboScan.py ${OUTDIR}/${filecounter}_dir/ -e fasta  -o ${OUTDIR}/${filecounter}_dir/scan/ -v 1

    python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ${OUTDIR}/${filecounter}_dir/scan/scannedScaffolds.gb  -o ${OUTDIR}/${filecounter}_dir/select/ -v 1

    python3.5 ~/GitHub/riboSeed/riboSeed/riboSnag.py ${OUTDIR}/${filecounter}_dir/scan/scannedScaffolds.gb  ${OUTDIR}/${filecounter}_dir/select/riboSelect_grouped_loci.txt  -o ${OUTDIR}/${filecounter}_dir/snag/ -v 1 -s ${i}
	cp ${OUTDIR}/${filecounter}_dir/snag/entropy_plot.pdf ${OUTDIR}/figures/${i}_entropy_plot.pdf
	cp ${OUTDIR}/${filecounter}_dir/snag/entropy_plot.png ${OUTDIR}/figures/${i}_entropy_plot.png

    filecounter=$((filecounter+1))
    IFS=$'\n\t'
done
