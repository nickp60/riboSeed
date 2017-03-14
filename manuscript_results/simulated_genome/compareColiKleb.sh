#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#

for i in "good" "bad";
do
    if [ $i == "good" ];
    then
	ref="NC_000913.3"
    else
	ref="CP003200.1"
    fi

# get the genome
echo "Downloading reference: ${ref}"
mkdir ${i}_ref
get_genomes.py -q $ref -o ./

# anotate the rDNAs

python3.5 ~/GitHub/riboSeed/riboSeed/riboScan.py ./${ref} .fasta -o ./${i}_ref/scan/
# cluster
python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ./${i}_ref/scan/scannedScaffolds.gb  -o ./${i}_ref/select/
# run riboSeed
python3.5 ~/GitHub/riboSeed/riboSeed/riboSeed.py -r ./${i}_ref/scan/scannedScaffolds.gb  -o ./${i}_ref/seed/ ./${i}_ref/select/riboSelect_grouped_loci.txt -F ./toyGenome/reads_1.fq -R ./toyGenome/reads_2.fq -z -v 1 -l 2000
done

# make copy of contigs renamed for mauve
echo "copying contigs to mauve dir"
cp ./good_ref/seed/final_de_fere_novo_assembly/contigs.fasta ./mauve/coli_de_fere_novo.fa

cp ./good_ref/seed/final_de_novo_assembly/contigs.fasta ./mauve/coli_de_novo.fa

cp ./bad_ref/seed/final_de_fere_novo_assembly/contigs.fasta ./mauve/kleb_de_fere_novo.fa
