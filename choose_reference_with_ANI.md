# Reference Selection
With any assembly that uses a reference, the choice of that reference is
crucial.  Here, we outline a protocol for using Average Nucleotide Identity via pyANI.

All the scripts are in the `scripts/ani_reference` dir.

Note: must have coreutils installed on OSX to provide the `gshuf` replacement for `shuf`.

## 1. Randomly Selecting Potential References
We use the list of prokaryote assemblies from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt). If you have it downloaded already, you can pass it as an arg to the get_n_genomes_script.

The following script will get 25 random E. coli genomes


```
/scripts/get_n_random_complete_genomes.sh -o "Escherichia coli" -n 25 > tmp_accessions
```
use the get_genomes.py script to download all of those into a directory called `potential_references`


```
mkdir potential_references
while read x; do echo "fetching $x" ; get_genomes.py -q $x -o potential_references ; done < tmp_accessions

```


## 2. Mini Assembly
Now, because we have a chicken/egg situation as we need an assembly to determine ANI to determine the best reference to make an assembly, we start with a minimal assembly with 1/10th of the reads using the `mini_assembly.sh` scsript.

```
./scripts/mini_assembly.sh ./tests/references/toy_reads1.fq ./tests/references/toy_reads2.fq TEST

# copy the resulting assembly scaffold into the same dir as the other genomes we selected
cp 2018-06-19_mini_assembly_TEST/spades/scaffolds.fasta ./potential_references/
```


## 3. Run pyANI
