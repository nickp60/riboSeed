#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#
FLANK=1000


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
if [ -e ${ref}.fasta ] ;
then
echo "using existing version of $ref"
else
get_genomes.py -q $ref -o ./
fi

# anotate the rDNAs

python3.5 ~/GitHub/riboSeed/riboSeed/riboScan.py ./${ref}.fasta -o ./${i}_ref/scan/
# cluster
python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ./${i}_ref/scan/scannedScaffolds.gb  -o ./${i}_ref/select/
# run riboSeed
python3.5 ~/GitHub/riboSeed/riboSeed/riboSeed.py -r ./${i}_ref/scan/scannedScaffolds.gb  -o ./${i}_ref/seed/ ./${i}_ref/select/riboSelect_grouped_loci.txt -F ./toyGenome/reads_1.fq -R ./toyGenome/reads_2.fq -i 3 -z -v 1 -l ${FLANK}  --keep_temps
done

# make copy of contigs renamed for mauve
echo "copying contigs to mauve dir"
cp ./toyGenome/coli_genome/scan/scannedScaffolds.gb ./mauve/reference.gb

cp ./good_ref/seed/final_de_fere_novo_assembly/contigs.fasta ./mauve/coli_de_fere_novo.fa

cp ./good_ref/seed/final_de_novo_assembly/contigs.fasta ./mauve/coli_de_novo.fa

cp ./bad_ref/seed/final_de_fere_novo_assembly/contigs.fasta ./mauve/kleb_de_fere_novo.fa
