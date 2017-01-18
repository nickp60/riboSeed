# RiboSeed Pipeline - especially quickstart
Sheesh, too impatient for the quickstart?

## Overview
No time, on to the examples!  *Schnell!*

## Minimal example
The pseudogenome was constructed from the 7 rDNAs separated by several kb of flanking DNA.  If can be found under ```./manuscript_results/simulated_genome/concatenated_seq.fasta```
Annotate rDNA regions with scanScaffolds.sh
```
riboScan.py ./manuscript_results/simulated_genome/ concatenated_seq.fasta -o ./scanned_output/ -k bac -t .5
```
Then, using riboSelect:

```
riboSelect.py ./scanned_output/scannedScaffolds.gb -o ./select_output/
```

No errors? riboSeed time!
```
riboSeed.py ./select_output/riboSelect_grouped_loci.txt \
            -F ./manuscript_results/simulated_genome/concatenated_seq1.fq \
            -R ./manuscript_results/simulated_genome/concatenated_seq2.fq \
            -r ./scanned_output/scannedScaffolds.gb \
            -o ./test_toy_reads/ -v 1  \
            --python2_7_exe /usr/bin/python2.7
```
