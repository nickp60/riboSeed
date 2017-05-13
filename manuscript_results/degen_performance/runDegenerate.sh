#!/bin/bash
set -euo pipefail
IFS=$'\n\t'


# process original
# only need to do this once;
# we arent doing indels so the coords are the same
mkdir ./degenerate_output/
mkdir ./degenerate_output/genomes/
mkdir ./degenerate_output/mauve/

python3.5 ~/GitHub/riboSeed/riboSeed/riboScan.py ./toyGenome/coli_genome/concatenated_seq.fasta -o ./degenerate_output/ref_scan/ -v 1


for i in 0 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1
do
    # intrduce the mutations
    python3.5  ~/GitHub/riboSeed/riboSeed/riboSim.py ./toyGenome/coli_genome/concatenated_seq.fasta -f ${i} -o ./degenerate_output/degenerate_${i}/ -v 1
    # copy the file
    cp ./degenerate_output/degenerate_${i}/concatenated_seq.fasta ./degenerate_output/genomes/
    # generate the gb file for the sequence
    python3.5 ~/GitHub/riboSeed/riboSeed/riboScan.py ./degenerate_output/degenerate_${i}/ -o ./degenerate_output/scan_${i}/ -e fasta -v 1
    # cluster
    python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ./degenerate_output/scan_${i}/scannedScaffolds.gb  -o ./degenerate_output/select_${i}/ -v 1
    head ./degenerate_output/select_${i}/riboSelect_grouped_loci.txt
    # run riboSeed
    python3.5 ~/GitHub/riboSeed/riboSeed/riboSeed.py -r ./degenerate_output/scan_${i}/scannedScaffolds.gb  -o ./degenerate_output/seed_${i}/ ./degenerate_output/select_${i}/riboSelect_grouped_loci.txt -F ./toyGenome/reads_1.fq -R ./toyGenome/reads_2.fq -i 3 -z -v 1 -l 1000 --keep_temps
    # move a copy of the annotated genome to the "mauve" directory for convenience
    cp ./degenerate_output/seed_${i}/final_de_fere_novo_assembly/contigs.fasta ./degenerate_output/mauve/de_fere_novo_${i}.fasta

done
