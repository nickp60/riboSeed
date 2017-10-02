# riboSeed - quickstart
Too Impatient for the quickstart? See our [Especially Quickstart Guide](./especiallyquickstart.md)

## Overview
In this guide, I will walk you through a typical experiment.

## Requirements

The following tools must be installed to complete this QuickStart guide

* [`sra-tools`/`sratoolkit`](https://github.com/ncbi/sra-tools)
* [`open_utils`](https://github.com/nickp60/open_utils)
* [`barrnap`](http://www.vicbioinformatics.com/software.barrnap.shtml)
* [`emboss`](http://www.ebi.ac.uk/Tools/emboss/)
* [`bwa`](http://bio-bwa.sourceforge.net/)
* [`riboSeed`](https://github.com/nickp60/riboSeed)

## Get the Data
Let's try using riboSeed on *Listeria monocytogenes*.  We will use the SRA toolkit's `fastq-dump` tool to get the data:

```
fastq-dump --split-files SRR5181495
```

That takes a while to download.

```
Read 618654 spots for SRR5181495
Written 618654 spots for SRR5181495
```

Let's use the *Listeria monocytogenes str. 4b F2365* genome as reference, NCBI accession [AE017262.2](https://www.ncbi.nlm.nih.gov/nuccore/AE017262).  Downloading this with the `get_genomes.py` tool installed along with [`pyutilsnrw` repository](https://github.com/nickp60/pyutilsnrw) (assuming this tool is in your `$PATH`):

```
$ get_genomes.py -q AE017262.2 -o ./
```

## `ribo run`
Now that we have the reference as a fasta file and our reads, lets run the pipeline:

```
$ ribo run ./AE017262.2.fasta -F ./SRR5181495_1.fastq -R ./SRR5181495_2.fastq  --cores 4 --threads 1 -v 1 --serialize  --keep_temps -o ./listeria_riboSeed_output/
```

## Evaluate
Use MAUVE to visualize the assemblies found in `./seed_output/final_de_novo_assembly` and `./seed_output/final_de_fere_novo_assembly`.  You should see that 4 out of the 6 regions were bridged.  The remaining region consists of a tandem rDNA repeat, which is a continuing problem for rDNA assembly.  Looking at the riboSelect output, we can see that it has chosen 6 clusters.  We can try to treat the tandem clusters as one (for a total of 5 clusters) by re-running `ribo run` with `--clusters 5`.
