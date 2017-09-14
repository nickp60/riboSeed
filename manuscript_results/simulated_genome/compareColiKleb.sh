#!/bin/bash
set -uo pipefail
IFS=$'\n\t'

#
FLANK=1000
if [ -z "$1" ]
then
    SEED=17
else
    SEED="$1"
fi
OUTDIR="./simulatedGenomeResults_${SEED}/"
mkdir ${OUTDIR}

##############################################################################
toyref="BA000007.2"
#ref="NC_000913.3"
# get the reference Sakai genome
if [ -e "${toyref}.fasta" ]
then
echo "using local copy of reference"
else
get_genomes.py -q ${toyref} -o ./
fi

# annotate regions
python3.5  ~/GitHub/riboSeed/riboSeed/riboScan.py ./${toyref}.fasta -o ${OUTDIR}toyGenome/scan/
# Cluster regions
python3.5  ~/GitHub/riboSeed/riboSeed/riboSelect.py ${OUTDIR}toyGenome/scan/scannedScaffolds.gb -o ${OUTDIR}toyGenome/select/
# extract regions with 5kb flanking
python3.5  ~/GitHub/riboSeed/riboSeed/riboSnag.py ${OUTDIR}toyGenome/scan/scannedScaffolds.gb  ${OUTDIR}toyGenome/select/riboSelect_grouped_loci.txt -o ${OUTDIR}toyGenome/snag/ -l 5000 --just_extract
# combine the extracted regions into a toy genome
python3.5 ~/GitHub/riboSeed/scripts/concatToyGenome.py ${OUTDIR}toyGenome/snag/ \*_riboSnag.fasta -o ${OUTDIR}toyGenome/coli_genome/
# generate reads from the toy genome simulating a MiSeq V3 run
~/bin/pirs-2.0.2/pirs simulate -m 300 -l 100 -x 30 -v 10 -o ${OUTDIR}/toyGenome/reads -B ~/bin/pirs-2.0.2/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz --compress -I ~/bin/pirs-2.0.2/Profiles/InDel_Profiles/phixv2.InDel.matrix -G ~/bin/pirs-2.0.2/Profiles/GC-depth_Profiles/humNew.gcdep_100.dat ${OUTDIR}/toyGenome/coli_genome/concatenated_seq.fasta -S ${SEED}


# annotate the toy genome for mauve visualization
python3.5  ~/GitHub/riboSeed/riboSeed/riboScan.py ${OUTDIR}toyGenome/coli_genome/concatenated_seq.fasta -o ${OUTDIR}toyGenome/coli_genome/scan/

##############################################################################

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
mkdir ${OUTDIR}${i}_ref
if [ -e ${ref}.fasta ] ;
then
echo "using existing version of $ref"
else
get_genomes.py -q $ref -o ./
fi

# anotate the rDNAs

python3.5 ~/GitHub/riboSeed/riboSeed/riboScan.py ./${ref}.fasta -o ${OUTDIR}${i}_ref/scan/
# cluster
python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ${OUTDIR}${i}_ref/scan/scannedScaffolds.gb  -o ${OUTDIR}${i}_ref/select/
# run riboSeed
python3.5 ~/GitHub/riboSeed/riboSeed/riboSeed.py -r ${OUTDIR}${i}_ref/scan/scannedScaffolds.gb  -o ${OUTDIR}${i}_ref/seed/ ${OUTDIR}${i}_ref/select/riboSelect_grouped_loci.txt -F ${OUTDIR}toyGenome/reads_100_300_1.fq.gz -R ${OUTDIR}toyGenome/reads_100_300_2.fq.gz -i 3 -c 2 -v 1 -l ${FLANK} --clean_temps
done

# make copy of contigs renamed for mauve
echo "copying contigs to mauve dir"
mkdir ${OUTDIR}mauve
cp ${OUTDIR}toyGenome/coli_genome/scan/scannedScaffolds.gb ${OUTDIR}mauve/reference.gb

cp ${OUTDIR}good_ref/seed/final_de_fere_novo_assembly/contigs.fasta ${OUTDIR}mauve/coli_de_fere_novo.fasta

cp ${OUTDIR}good_ref/seed/final_de_novo_assembly/contigs.fasta ${OUTDIR}mauve/coli_de_novo.fasta

cp ${OUTDIR}bad_ref/seed/final_de_fere_novo_assembly/contigs.fasta ${OUTDIR}mauve/kleb_de_fere_novo.fasta

python3.5 ~/GitHub/riboSeed/riboSeed/riboSketch.py ${OUTDIR}mauve/ -o ${OUTDIR}riboSketch/  --names "Artificial Genome, De fere novo, De novo, Kleb. de fere novo"

python3.5 ~/GitHub/riboSeed/riboSeed/riboScore.py ${OUTDIR}mauve -o ${OUTDIR}riboScore/
