#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

ref="BA000007.2"
#ref="NC_000913.3"
# get the reference Sakai genome
if [ -e "${ref}.fasta" ]
then
echo "using local copy of reference"
else
get_genomes.py -q ${ref} -o ./
fi

# annotate regions
python3.5  ~/GitHub/riboSeed/riboSeed/riboScan.py ./${ref}.fasta -o ./toyGenome/scan/
# Cluster regions
python3.5  ~/GitHub/riboSeed/riboSeed/riboSelect.py ./toyGenome/scan/scannedScaffolds.gb -o ./toyGenome/select/
# extract regions with 5kb flanking
python3.5  ~/GitHub/riboSeed/riboSeed/riboSnag.py ./toyGenome/scan/scannedScaffolds.gb  ./toyGenome/select/riboSelect_grouped_loci.txt -o ./toyGenome/snag/ -l 5000 --just_extract
# combine the extracted regions into a toy genome
python3.5 ~/GitHub/riboSeed/scripts/concatToyGenome.py ./toyGenome/snag/ \*_riboSnag.fasta -o ./toyGenome/coli_genome/
# generate reads from the toy genome simulating a MiSeq V3 rub

# ~/bin/art_bin_MountRainier/art_illumina -ss HS20 -i ./toyGenome/coli_genome/concatenated_seq.fasta -p -l 100 -f 30 -m 300 -s 10 -o ./toyGenome/reads_100_300_ -rs 3
# gzip -c ./toyGenome/reads_100_300_1.fq > ./toyGenome/reads_100_300_1.fq.gz
# gzip -c ./toyGenome/reads_100_300_2.fq > ./toyGenome/reads_100_300_2.fq.gz


~/bin/pIRS_111/pirs simulate -i ./toyGenome/coli_genome/concatenated_seq.fasta -m 300 -l 100 -x 30 -v 10 -o ./toyGenome/reads


# annotate the toy genome for mauve visualization
python3.5  ~/GitHub/riboSeed/riboSeed/riboScan.py ./toyGenome/coli_genome/concatenated_seq.fasta -o ./toyGenome/coli_genome/scan/
