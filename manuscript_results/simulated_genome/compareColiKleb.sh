#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#
FLANK=1000


for i in "good" "bad";
do
    if [ $i == "good" ];
    then
	# ref="BA000007.2"
	ref="NC_000913.3"
	# ref="NZ_CP008957.1"
	# ref="NC_011751.1"
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
python3.5 ~/GitHub/riboSeed/riboSeed/riboSeed.py -r ./${i}_ref/scan/scannedScaffolds.gb  -o ./${i}_ref/seed/ ./${i}_ref/select/riboSelect_grouped_loci.txt -F ./toyGenome/reads_100_300_1.fq.gz -R ./toyGenome/reads_100_300_2.fq.gz -i 3 -z -v 1 -l ${FLANK}
done

# make copy of contigs renamed for mauve
echo "copying contigs to mauve dir"
mkdir mauve
cp ./toyGenome/coli_genome/scan/scannedScaffolds.gb ./mauve/reference.gb

cp ./good_ref/seed/final_de_fere_novo_assembly/contigs.fasta ./mauve/coli_de_fere_novo.fa

cp ./good_ref/seed/final_de_novo_assembly/contigs.fasta ./mauve/coli_de_novo.fa

cp ./bad_ref/seed/final_de_fere_novo_assembly/contigs.fasta ./mauve/kleb_de_fere_novo.fa


python3.5 ~/GitHub/riboSeed/scripts/plotMauveBetter.py ./mauve/reference.gb ./mauve/coli_de_fere_novo.fa ./mauve/coli_de_novo.fa ./mauve/kleb_de_fere_novo.fa -o ./ordered_plots/  --names "Artificial Genome, De fere novo, De novo, Kleb. de fere novo"
# java -Xmx500m -cp ~/mauve_snapshot_2015-02-13/Mauve.jar  org.gel.mauve.contigs.ContigOrderer -output de_fere -ref ./mauve/reference.gb -draft ./mauve/coli_de_fere_novo.fa

# java -Xmx500m -cp ~/mauve_snapshot_2015-02-13/Mauve.jar  org.gel.mauve.contigs.ContigOrderer -output de_novo -ref ./mauve/reference.gb -draft ./mauve/coli_de_novo.fa

# java -Xmx500m -cp ~/mauve_snapshot_2015-02-13/Mauve.jar  org.gel.mauve.contigs.ContigOrderer -output de_fere_kleb -ref ./mauve/reference.gb -draft ./mauve/kleb_de_fere_novo.fa
