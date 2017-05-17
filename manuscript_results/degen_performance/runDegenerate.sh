#!/bin/bash
# version 0.2.0
IFS=$'\n\t'


# process original
# only need to do this once;
# we arent doing indels so the coords are the same
SEED=2
OUTDIRBASE="./"
OUTDIR="${OUTDIRBASE}`date +%F`_degenerate_output_${SEED}"


mkdir ${OUTDIR}
mkdir ${OUTDIR}/genomes/
mkdir ${OUTDIR}/mauve/

LOGFILE=${OUTDIR}/log.txt


python3.5 ~/GitHub/riboSeed/riboSeed/riboScan.py ./toyGenome/coli_genome/ -o ${OUTDIR}/ref_scan/ -v 1 -e fasta &>> ${LOGFILE}

python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ${OUTDIR}/ref_scan/scannedScaffolds.gb  -o ${OUTDIR}/ref_select/ -v 1  &>> ${LOGFILE}

for i in 0.0001 0.00025 0.0005 0.00075 0.001 0.0025 0.005 0.0075 0.01 0.05
do
    echo "Processing mutation frequency of ${i}"
    # intrduce the mutations
    python3.5  ~/GitHub/riboSeed/riboSeed/riboSim.py ./toyGenome/coli_genome/concatenated_seq.fasta -f ${i} -o ${OUTDIR}/degenerate_${i}/ -v 1 --seed ${SEED} &>> ${LOGFILE}
    # copy the file
    cp ${OUTDIR}/degenerate_${i}/concatenated_seq.fasta ${OUTDIR}/genomes/${i}_reference.fasta
    # generate the gb file for the sequence
    python3.5 ~/GitHub/riboSeed/riboSeed/riboScan.py ${OUTDIR}/degenerate_${i}/ -o ${OUTDIR}/scan_${i}/ -e fasta -v 1  &>> ${LOGFILE}
    # cluster
    # python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ${OUTDIR}/scan_${i}/scannedScaffolds.gb  -o ${OUTDIR}/select_${i}/ -v 1
    # head ${OUTDIR}/select_${i}/riboSelect_grouped_loci.txt
    # run riboSeed
    python3.5 ~/GitHub/riboSeed/riboSeed/riboSeed.py -r ${OUTDIR}/scan_${i}/scannedScaffolds.gb  -o ${OUTDIR}/seed_${i}/ ${OUTDIR}/ref_select/riboSelect_grouped_loci.txt -F ./toyGenome/reads_1.fq -R ./toyGenome/reads_2.fq -i 3  -v 1 -l 1000 --keep_temps  &>> ${LOGFILE}
    # move a copy of the annotated genome to the "mauve" directory for convenience
    cp ${OUTDIR}/seed_${i}/final_de_fere_novo_assembly/contigs.fasta ${OUTDIR}/mauve/de_fere_novo_${i}.fasta

done
