# riboSeed - especially quickstart
Sheesh, too impatient for the quickstart?

## Overview
No time, on to the examples!  *Schnell!*

## Minimal example
The pseudogenome was constructed from the 7 rDNAs separated by several kb of flanking DNA.  If can be found under `./riboSeed/integration_test/concatenated_seq.fasta` in the git repo.  If you have installed using setuptools, the `integration_test` folder will be installed in the site-packages dir, such as `/venv-riboSeed/lib/python3.5/site-packages/riboSeed/integration_data/`.

Two read files can be found in the same directory.

Run the analysis using the built in `ribo try` command.
```
ribo test

# This is running the following command under the hood,
# ribo run -r ./riboSeed/integration_data/concatenated_seq.fasta -F ./riboSeed/integration_data/test_reads1.fq -R ./riboSeed/integration_data/test_reads2.fq -o ./test1/ -v 1
```
