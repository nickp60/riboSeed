# riboSeed - quickstart
Too Impatient for the quickstart? See our [Especially Quickstart Guide](./especiallyquickstart.md)

## Overview
In this guide, I will walk you through a typical experiment.

## Requirements

```
conda install riboseed
```

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

Let's use the *Listeria monocytogenes str. 4b F2365* genome as reference, NCBI accession [AE017262.2](https://www.ncbi.nlm.nih.gov/nuccore/AE017262).  Download this :
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/285/GCF_000008285.1_ASM828v1/GCF_000008285.1_ASM828v1_genomic.fna.gz
# unzip it
gunzip GCF_000008285.1_ASM828v1_genomic.fna.gz
# rename it something slightly shorter
mv GCF_000008285.1_ASM828v1_genomic.fna  AE017262.2.fasta
```

## `ribo run`
Now that we have the reference as a fasta file and our reads, lets run the pipeline:

```
$ ribo run -r ./AE017262.2.fasta -F ./SRR5181495_1.fastq -R ./SRR5181495_2.fastq  --cores 4 --threads 1 -v 1 --serialize -o ./listeria_riboSeed_output/
```

## Evaluate
Use MAUVE to visualize the assemblies found in `./seed_output/final_de_novo_assembly` and `./seed_output/final_de_fere_novo_assembly`.  You should see that 4 out of the 6 regions were bridged.  The remaining region consists of a tandem rDNA repeat, which is a continuing problem for rDNA assembly.  Looking at the riboSelect output, we can see that it has chosen 6 clusters.  We can try to treat the tandem clusters as one (for a total of 5 clusters) by re-running `ribo run` with `--clusters 5`.
