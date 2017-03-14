#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# get the reference Sakai genome
get_genomes.py -q BA000007.2 -o ./

# annotate regions
python3.5  ~/GitHub/riboSeed/riboSeed/riboScan.py ./BA000007.2 .fasta -o ./toyGenome/scan/
# Cluster regions
python3.5  ~/GitHub/riboSeed/riboSeed/riboSelect.py ./toyGenome/scan/scannedScaffolds.gb -o ./toyGenome/select/
# extract regions with 5kb flanking
python3.5  ~/GitHub/riboSeed/riboSeed/riboSnag.py ./toyGenome/scan/scannedScaffolds.gb  ./toyGenome/select/riboSelect_grouped_loci.txt -o ./toyGenome/snag/ -l 5000
# combine the extracted regions into a toy genome
python3.5 ~/GitHub/riboSeed/scripts/concatToyGenome.py ./toyGenome/snag/ \*_riboSnag.fasta -o ./toyGenome/coli_genome/
# generate reads from the toy genome simulating a MiSeq V3 rub
~/bin/art_bin_MountRainier/art_illumina -ss MSv3 -i ./toyGenome/coli_genome/concatenated_seq.fasta -p -l 250 -f 20 -m 300 -s 10 -o ./toyGenome/reads_
# annotate the toy genome for mauve visualization
python3.5  ~/GitHub/riboSeed/riboSeed/riboScan.py ./toyGenome/coli_genome/concatenat .fasta -o ./toyGenome/coli_genome/scan/


# move a copy of the annotated genome to the "mauve" directory for convenience
mkdir mauve
cp ./toyGenome/coli_genome/scan/scannedScaffolds.gb ./mauve/concatenated_coli.gb
