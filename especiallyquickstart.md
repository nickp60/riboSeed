# RiboSeed Pipeline - especially quickstart
Sheesh, too impatient for the quickstart?

## Overview
No time, on to the examples!  *Schnell!*
## Minimal example
The pseudogenome was constructed from the 7 rDNAs separated by several kb of flanking DNA.  If can be found under ```./sample_data/pseudo/combined_colapsed.fa```
Annotate rDNA regions with scanScaffolds.sh
```
 ~/GitHub/riboSeed/scripts/scanScaffolds.sh ./sample_data/toy_set/pseudata/ fa ./scan_toy/ bac .5
```
Then, using riboSelect:

```
python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ./scan_toy/scannedScaffolds.gb -o ./test_output/
```

No errors? riboSeed time!
```
python3.5 /home/nicholas/GitHub/riboSeed/riboSeed/riboSeed2.py ./test_select/riboSelect_grouped_loci.txt -F ./sample_data/toy_set/toy_reads1.fq -R ./sample_data/toy_set/toy_reads2.fq -r ./scan/scannedScaffolds.gb -o ./test_toy_reads/ -v 1 -i 1 --min_assembly_len 6000  --python2_7_exe /usr/bin/python2.7 -t 1
```
