# RiboSeed Pipeline - especially quickstart
Sheesh, too impatient for the quickstart?

## Overview
No time, on to the examples!  *Schnell!*

## Minimal example
The pseudogenome was constructed from the 7 rDNAs separated by several kb of flanking DNA.  If can be found under ```./manuscript_results/simulated_genome/concatenated_seq.fasta```
Annotate rDNA regions with scanScaffolds.sh
```
./scripts/scanScaffolds.sh sh ./sample_data/toy_set/pseudata/ seq.fa ./scanned_output/ bac .5
```
Then, using riboSelect:

```
python3.5 ./riboSeed/riboSelect.py ./scanned_output/scannedScaffolds.gb -o ./select_output/
```

No errors? riboSeed time!
```
python3.5 ./riboSeed/riboSeed.py ./select_output/riboSelect_grouped_loci.txt -F ./manuscript_results/simulated_genome/concatenated_seq1.fq -R ./manuscript_results/simulated_genome/concatenated_seq2.fq -r ./scanned_output/scannedScaffolds.gb -o ./test_toy_reads/ -v 1  --python2_7_exe /usr/bin/python2.7
```
