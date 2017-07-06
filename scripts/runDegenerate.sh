#!/bin/bash
set -euo pipefail
# version 0.2.0
IFS=$'\n\t'


# process original
# only need to do this once;
# we arent doing indels so the coords are the same

if [ -z "$1" ]
then
    SEED=17
else
    SEED="$1"
fi

OUTDIRBASE="./"
OUTDIR="${OUTDIRBASE}`date +%F`_degenerate_output_${SEED}"

mkdir ${OUTDIR}
mkdir ${OUTDIR}/toyGenome/
mkdir ${OUTDIR}/genomes/
mkdir ${OUTDIR}/mauve/

LOGFILE=${OUTDIR}/log.txt

echo "Running degenerate pipeline with the following settings"
echo "OUTDIR: ${OUTDIR}"
echo "SEED: ${SEED}"

# get the reference Sakai genome
if [ -e BA000007.2.fasta ]
then
echo "using local copy of reference"
else
get_genomes.py -q BA000007.2 -o ./
fi


echo " annotate regions in full genome"
python3.5  ~/GitHub/riboSeed/riboSeed/riboScan.py ./BA000007.2.fasta -o ${OUTDIR}/toyGenome/scan/ &>> ${LOGFILE}
echo " Cluster regions"
python3.5  ~/GitHub/riboSeed/riboSeed/riboSelect.py ${OUTDIR}/toyGenome/scan/scannedScaffolds.gb -o ${OUTDIR}/toyGenome/select/ &>> ${LOGFILE}
echo "extract regions with 5kb flanking"
python3.5  ~/GitHub/riboSeed/riboSeed/riboSnag.py ${OUTDIR}/toyGenome/scan/scannedScaffolds.gb  ${OUTDIR}/toyGenome/select/riboSelect_grouped_loci.txt -o ${OUTDIR}/toyGenome/snag/ -l 5000 -v 1 --just_extract &>> ${LOGFILE}
echo "combine the extracted regions into a toy genome"
python3.5 ~/GitHub/riboSeed/scripts/concatToyGenome.py ${OUTDIR}/toyGenome/snag/ \*_riboSnag.fasta -o ${OUTDIR}/toyGenome/coli_genome/ &>> ${LOGFILE}
echo " generate reads from the toy genome"
# ~/bin/art_bin_MountRainier/art_illumina -ss MSv3 -i ${OUTDIR}/toyGenome/coli_genome/concatenated_seq.fasta -p -l 250 -f 20 -m 300 -s 10 -rs ${SEED} -o ${OUTDIR}/toyGenome/reads_ &>> ${LOGFILE}

# ~/bin/pIRS_111/pirs simulate -i ${OUTDIR}/toyGenome/coli_genome/concatenated_seq.fasta -m 300 -l 100 -x 30 -v 10 -o ${OUTDIR}/toyGenome/reads &>> ${LOGFILE}
~/bin/pirs-2.0.2/pirs simulate -m 300 -l 100 -x 30 -v 10 -o ${OUTDIR}/toyGenome/reads -B ~/bin/pirs-2.0.2/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz --compress -I ~/bin/pirs-2.0.2/Profiles/InDel_Profiles/phixv2.InDel.matrix -G ~/bin/pirs-2.0.2/Profiles/GC-depth_Profiles/humNew.gcdep_100.dat ${OUTDIR}/toyGenome/coli_genome/concatenated_seq.fasta &>> ${LOGFILE}


# annotate the toy genome for mauve visualization
python3.5  ~/GitHub/riboSeed/riboSeed/riboScan.py ${OUTDIR}/toyGenome/coli_genome/ -o  ${OUTDIR}/ref_scan/ -v 1 -e fasta &>> ${LOGFILE}



# move a copy of the annotated genome to the "mauve" directory for convenience
cp ${OUTDIR}/ref_scan/scannedScaffolds.gb ${OUTDIR}/mauve/reference.gb


# cluster regions in toy genome
python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ${OUTDIR}/ref_scan/scannedScaffolds.gb  -o ${OUTDIR}/ref_select/ -v 1  &>> ${LOGFILE}

SUBSEEDS=(`~/GitHub/riboSeed/scripts/seedRand.py $SEED 10`)
FREQS=( 0.0 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 )
for ((i=0; i<${#SUBSEEDS[@]}; ++i));
do

    SUBSEED=${SUBSEEDS[i]}
    FREQ=${FREQS[i]}
# do
#     echo $i
#     echo "${SUBSEEDS[i]}"
#     echo "${FREQS[i]}"
# done
# exit 0
# for i in 0.0 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1
# do
    echo "Processing mutation frequency of ${FREQ}"
    ###################################################################
    # make a new reference
    mkdir ${OUTDIR}/mutated_regions_${FREQ}/
    for region in ${OUTDIR}/toyGenome/snag/*_riboSnag.fasta
    do
	echo "Mutating ${region}"
	# intrduce the mutations
	fname_ext=${region##*/}
	fname=${fname_ext%.*}
	python3.5  ~/GitHub/riboSeed/riboSeed/riboSim.py ${region} -e 5000  -f ${FREQ} -o ${OUTDIR}/degenerate_${FREQ}_${fname}/ -v 1 --seed ${SUBSEED} &>> ${LOGFILE}
	# copy results to a directory
	cp ${OUTDIR}/degenerate_${FREQ}_${fname}/${fname_ext} ${OUTDIR}/mutated_regions_${FREQ}/
    done
    # combine the mutated extracted regions into a toy genome
    python3.5 ~/GitHub/riboSeed/scripts/concatToyGenome.py ${OUTDIR}/mutated_regions_${FREQ}/ \*_riboSnag.fasta -o ${OUTDIR}/toyGenome/mutated_genome_${FREQ}/ &>> ${LOGFILE}
    ###################################################################
    # copy the file
    cp ${OUTDIR}/toyGenome/mutated_genome_${FREQ}/concatenated_seq.fasta ${OUTDIR}/genomes/${FREQ}_reference.fasta
    # generate the gb file for the sequence
    python3.5 ~/GitHub/riboSeed/riboSeed/riboScan.py ${OUTDIR}/toyGenome/mutated_genome_${FREQ}/ -o ${OUTDIR}/scan_${FREQ}/ -e fasta -v 1  &>> ${LOGFILE}
    # run riboSeed
    python3.5 ~/GitHub/riboSeed/riboSeed/riboSeed.py -r ${OUTDIR}/scan_${FREQ}/scannedScaffolds.gb  -o ${OUTDIR}/seed_${FREQ}/ ${OUTDIR}/ref_select/riboSelect_grouped_loci.txt -F ${OUTDIR}/toyGenome/reads_100_300_1.fq.gz -R ${OUTDIR}/toyGenome/reads_100_300_2.fq.gz -i 3  -v 1 -l 1000 --clean_temps  &>> ${LOGFILE}
    # move a copy of the annotated genome to the "mauve" directory for convenience
    cp ${OUTDIR}/seed_${FREQ}/final_de_fere_novo_assembly/contigs.fasta ${OUTDIR}/mauve/de_fere_novo_${FREQ}.fasta

    python3.5 ~/GitHub/riboSeed/riboSeed/riboScore.py ${OUTDIR}/seed_${FREQ}/mauve/

done
