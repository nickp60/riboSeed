#!/bin/bash
for len in 50 75 100 150 200
# for each coverage level, generate some MiSeq V3 reads
do
for cov in 10 20 50 100;
do
    echo "trying to get len ${len} and cov ${cov}"
    scp ppserver:/mnt/shared/scratch/nw42839/simulated_reads_coli/2017-05-18_simulation/sim_len_${len}_cov_${cov}/_seed/combined_quast_report.tsv ./sim_len_${len}_cov_${cov}.csv

done
done
