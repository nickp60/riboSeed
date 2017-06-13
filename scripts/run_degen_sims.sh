#!/bin/bash
mkdir mauve_results
DATE=`date +%F`
for i in 9 13 17 23 27 31 42 49
do
 ../scripts/runDegenerate.sh ${i}
done

for i in 9 13 17 23 27 31 42 49
do
    progressiveMauve --output=mauve_results/seed_${i}.xfma ./${DATE}_degenerate_output_${i}/mauve/de_fere_novo_0.02.fasta ./${DATE}_degenerate_output_${i}/mauve/de_fere_novo_0.03.fasta ./${DATE}_degenerate_output_${i}/mauve/de_fere_novo_0.04.fasta ./${DATE}_degenerate_output_${i}/mauve/de_fere_novo_0.05.fasta ./${DATE}_degenerate_output_${i}/mauve/de_fere_novo_0.06.fasta ./${DATE}_degenerate_output_${i}/mauve/de_fere_novo_0.07.fasta ./${DATE}_degenerate_output_${i}/mauve/de_fere_novo_0.08.fasta ./${DATE}_degenerate_output_${i}/mauve/de_fere_novo_0.09.fasta ./${DATE}_degenerate_output_${i}/mauve/de_fere_novo_0.0.fasta ./${DATE}_degenerate_output_${i}/mauve/de_fere_novo_0.1.fasta ./${DATE}_degenerate_output_${i}/mauve/reference.gb



done
