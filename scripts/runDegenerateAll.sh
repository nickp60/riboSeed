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
OUTDIR="${OUTDIRBASE}`date +%F`_degenerate_output_ALL_${SEED}"

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
echo "running:"
~/bin/pirs-2.0.2/pirs simulate -m 300 -l 100 -x 30 -v 10 -o ${OUTDIR}/toyGenome/reads -B ~/bin/pirs-2.0.2/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz --compress -I ~/bin/pirs-2.0.2/Profiles/InDel_Profiles/phixv2.InDel.matrix -G ~/bin/pirs-2.0.2/Profiles/GC-depth_Profiles/humNew.gcdep_100.dat ${OUTDIR}/toyGenome/coli_genome/concatenated_seq.fasta &>> ${LOGFILE}


# annotate the toy genome for mauve visualization
python3.5  ~/GitHub/riboSeed/riboSeed/riboScan.py ${OUTDIR}/toyGenome/coli_genome/ -o  ${OUTDIR}/ref_scan/ -v 1 -e fasta &>> ${LOGFILE}



# move a copy of the annotated genome to the "mauve" directory for convenience
cp ${OUTDIR}/ref_scan/scannedScaffolds.gb ${OUTDIR}/mauve/reference.gb


# cluster regions in toy genome
python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ${OUTDIR}/ref_scan/scannedScaffolds.gb  -o ${OUTDIR}/ref_select/ -v 1  &>> ${LOGFILE}

for i in 0.0 0.0001 0.00025 0.0005 0.00075 0.001 0.0025 0.005 0.0075 0.01 0.05
	 # for i in  0.00075
	 # for i in 0.0 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1
do
    echo "Processing mutation frequency of ${i}"
    ###################################################################
    # make a new reference
    mkdir ${OUTDIR}/mutated_regions_${i}/
    for region in ${OUTDIR}/toyGenome/snag/*_riboSnag.fasta
    do
	echo "Mutating ${region}"
	# intrduce the mutations
	fname_ext=${region##*/}
	fname=${fname_ext%.*}
	# echo $fname
	# echo $fname_ext
	python3.5  ~/GitHub/riboSeed/riboSeed/riboSim.py ${region}  -f ${i} -o ${OUTDIR}/degenerate_${i}_${fname}/ -v 1 --seed ${SEED} &>> ${LOGFILE}
	# copy results to a directory
	cp ${OUTDIR}/degenerate_${i}_${fname}/${fname_ext} ${OUTDIR}/mutated_regions_${i}/
    done
    # combine the mutated extracted regions into a toy genome
    python3.5 ~/GitHub/riboSeed/scripts/concatToyGenome.py ${OUTDIR}/mutated_regions_${i}/ \*_riboSnag.fasta -o ${OUTDIR}/toyGenome/mutated_genome_${i}/ &>> ${LOGFILE}
    ###################################################################
    # copy the file
    cp ${OUTDIR}/toyGenome/mutated_genome_${i}/concatenated_seq.fasta ${OUTDIR}/genomes/${i}_reference.fasta
    # generate the gb file for the sequence
    python3.5 ~/GitHub/riboSeed/riboSeed/riboScan.py ${OUTDIR}/toyGenome/mutated_genome_${i}/ -o ${OUTDIR}/scan_${i}/ -e fasta -v 1  &>> ${LOGFILE}
    # cluster
    # python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ${OUTDIR}/scan_${i}/scannedScaffolds.gb  -o ${OUTDIR}/select_${i}/ -v 1
    # head ${OUTDIR}/select_${i}/riboSelect_grouped_loci.txt
    # run riboSeed
    python3.5 ~/GitHub/riboSeed/riboSeed/riboSeed.py -r ${OUTDIR}/scan_${i}/scannedScaffolds.gb  -o ${OUTDIR}/seed_${i}/ ${OUTDIR}/ref_select/riboSelect_grouped_loci.txt -F ${OUTDIR}/toyGenome/reads_100_300_1.fq.gz -R ${OUTDIR}/toyGenome/reads_100_300_2.fq.gz -i 3  -v 1 -l 1000 --clean_temps  &>> ${LOGFILE}
    # move a copy of the annotated genome to the "mauve" directory for convenience


    if [ -e ${OUTDIR}/seed_${i}/final_de_fere_novo_assembly/contigs.fasta ]
    then
	cp ${OUTDIR}/seed_${i}/final_de_fere_novo_assembly/contigs.fasta ${OUTDIR}/mauve/de_fere_novo_${i}.fasta
    fi
    
    if [ -e ${OUTDIR}/seed_${i}/final_de_novo_assembly/contigs.fasta ]
    then
	cp ${OUTDIR}/seed_${i}/final_de_novo_assembly/contigs.fasta ${OUTDIR}/mauve/de_novo_${i}.fasta
    fi
    
done
