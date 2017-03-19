# RiboSeed Pipeline - quickstart
Too Impatient for the quickstart? See our [Especially Quickstart Guide](./especiallyquickstart.md)

## Overview
In this guide, I will walk you through a typical experiment.

## Requirements

The following tools must be installed to complete this QuickStart guide

(On a mac? Try the included `install_deps.sh` script to download all the required tools and add them to your PATH)


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

## riboScan
Now that we have the reference as a fasta file, lets use `riboScan.py` (essentially barrnap) to annotate the rDNAs in the genome:

```
$ riboScan ./AE017262.2.fasta -o ./scan/
```

## riboSelect

Now that we have our genbank file, lets try to run riboSelect (assuming that `riboSelect.py` is in your `$PATH`) to pull out the regions of interest.

```
$ riboSelect.py ./scan/scannedScaffolds.gb -o ./select_output/
2017-01-18 09:31:22 - INFO - Initializing logger
2017-01-18 09:31:22 - INFO - Current usage:
./AE017262.2.gb -o ./scanned_output/

Loaded 1 records
2017-01-18 09:31:22 - INFO - Processing AE017262.2

2017-01-18 09:31:22 - INFO - hits in AE017262.2

#$ FEATURE rRNA
# Generated cluters for 0 on 20170118
AE017262.2 LMOf2365_2904:LMOf2365_2905:LMOf2365_2906
AE017262.2 LMOf2365_2927:LMOf2365_2928:LMOf2365_2931
AE017262.2 LMOf2365_2851:LMOf2365_2861:LMOf2365_2862:LMOf2365_2863
AE017262.2 LMOf2365_2936:LMOf2365_2937:LMOf2365_2940
AE017262.2 LMOf2365_2849:LMOf2365_2850
AE017262.2 LMOf2365_2900:LMOf2365_2901:LMOf2365_2902

```

Now we have our cluster file, reference gb, and data, we can apply `riboSeed.py`...

## riboSeed

```
riboSeed.py ./select_output/riboSelect_grouped_loci.txt -o ./seed_output/ -F ./SRR5181495_1.fastq -R ./SRR5181495_2.fastq -r ./AE017262.2.gb -c 4 -t 1 -v 1 --serialize  --keep_temps
```

## Evaluate
Use MAUVE to visualize the assemblies found in `./seed_output/final_de_novo_assembly` and `./seed_output/final_de_fere_novo_assembly`.  You should see that 4 out of the 6 regions were bridged.  The remaining region consists of a tandem rDNA repeat, which is a continuing problem for rDNA assembly.  Looking at the riboSelect output, we can see that it has chosen 6 clusters.  We can try to treat the tandem clusters as one (for a total of 5 clusters) by re-running riboSelect with the `-c 5`.
