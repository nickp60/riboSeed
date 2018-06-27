# Reference Selection via ANI
With any assembly that uses a reference, the choice of that reference is
crucial. Here, we outline a protocol for using Average Nucleotide Identity via pyANI.

All the subscripts are in the `scripts/select_ref_by_ANI` dir, and the main runner script can be found at `scripts/run_ANI.sh`

# Required Tools
## `shuf`
You must have coreutils installed on OSX to provide the `gshuf` replacement for `shuf`.

## `pyani`
`pyani` makes it easy to perform average nucleotide identity computations.  install with either `conda install pyani` (prefered), or with `pip install pyani`.

# `run_ANI.sh`: Helper Script Usage
This script helps automate this process.  It takes 5 mandatory args:


- `-e` experiement name
- `-f` Forward Reads
- `-r` Reverse Reads
- `-n` number of strains
- `-o` organism name

 And 2 optional arguments.  This allows you to reuse the same prokaryotes.txt file from NCBI and the same comparison genomes (rather than having to download them fresh each time)

- `-p` prokaryotes.txt file
- `-g` directory of genomes of interest

This returns the string of the file name of the best reference, as well as writes that to a file in the output directory called `best_reference`.


What's it doing under the hood?

# Description of Components
## 1. Randomly Selecting Potential References
We use the list of prokaryote assemblies from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt). If you have it downloaded already, you can pass it as an arg to the get_n_genomes_script.

The following script will get 25 random E. coli genomes


```
./scripts/select_ref_by_ANI/get_n_random_complete_genomes.sh -o "Escherichia coli" -n 25 > tmp_accessions
```
use the get_genomes.py script to download all of those into a directory called `potential_references`


```
mkdir potential_references
while read x
do
  echo "fetching $x"
  get_genomes.py -q $x -o potential_references
done < tmp_accessions

```


## 2. Mini Assembly
Now, because we have a chicken/egg situation as we need an assembly to determine ANI to determine the best reference to make an assembly, we start with a minimal assembly with 1/10th of the reads using the `mini_assembly.sh` scsript.

```
./scripts/mini_assembly.sh ./tests/references/toy_reads1.fq ./tests/references/toy_reads2.fq TEST

# copy the resulting assembly scaffold into the same dir as the other genomes we selected
cp 2018-06-19_mini_assembly_TEST/spades/scaffolds.fasta ./potential_references/
```


## 3. Run pyANI

```
average_nucleotide_identity.py -i potential_references -g -o ./pyani

```
Look at the resulting `pyani/ANIm_percentage_identity.tab` file; the best hit will be the one in the column/row for "contigs" with the closest score to 1.
