# RiboSeed Pipeline - quickstart
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
Let's try using riboSeed on S. aureus SS_0588, from https://www.ncbi.nlm.nih.gov/bioproject/PRJNA312437.  We will use the SRA toolkit's `fastq-dunp` tool to get the data:

```
fastq-dump --split-files SRR4360364
```

That takes a while to download.

```
Read 2231775 spots for SRR4360364
Written 2231775 spots for SRR4360364
```

Let's use the MRSA252 reference, NCBI BX571856.1.  Downloading this with the `get_genomes.py` tool installed along with [`pyutilsnrw` repository](https://github.com/nickp60/pyutilsnrw) (assuming this tool is in your `$PATH`):

```
$ ../open_utils/get_genomes/get_genomes.py -q BX571856.1 -o ./ -f gb
```


## RiboSelect

Now that we have our genbank file, lets try to run riboSelect (assuming that `riboSelect.py` is in your `$PATH`) to pull out the regions of interest.

Because this genome is pretty old (2004), we will have to reannotate the rDNA loci.  To do this, we will download the fasta version of the genome, and use the `riboScan.py` script that wraps `barrnap` and `seqret` (we assume that this script is in your `$PATH`).

```
$ get_genomes.py -q BX571856.1 -o ./ -f fasta
$ riboScan.py ./ .fasta -o ./scanned_output/ -k bac
example: $ barrnap
using default threshold of .5 identity
Processing BX571856.1.fasta, item 1 out of 1
[12:07:46] This is barrnap 0.7
[12:07:46] Written by Torsten Seemann <torsten.seemann@gmail.com>
[12:07:46] Obtained from https://github.com/Victorian-Bioinformatics-Consortium/barrnap
[12:07:46] Detected operating system: linux
[12:07:46] Using HMMER binary: /home/nicholas/barrnap/bin/../binaries/linux/nhmmer
[12:07:46] Will use 8 threads
[12:07:46] Setting evalue cutoff to 1e-06
[12:07:46] Will tag genes  < 0.8 of expected length.
[12:07:46] Will reject genes < .5 of expected length.
[12:07:46] Using database: /home/nicholas/barrnap/bin/../db/bac.hmm
[12:07:46] Scanning ./scanned_output/BX571856.1.fasta_barrnap_renamed.fasta for bac rRNA genes... please wait
[12:07:46] Command: /home/nicholas/barrnap/bin/../binaries/linux/nhmmer --cpu 8 -E 1e-06 --w_length 3878  -o /dev/null --tblout /dev/stdout \/home\/nicholas\/barrnap\/bin\/\.\.\/db\/bac\.hmm \.\/scanned_output\/BX571856\.1\.fasta_barrnap_renamed\.fasta
[12:07:49] Found: 16S_rRNA BX571856.1 L=1550/1585 514254..515803 + 16S ribosomal RNA
[12:07:49] Found: 16S_rRNA BX571856.1 L=1550/1585 559981..561530 + 16S ribosomal RNA
[12:07:49] Found: 16S_rRNA BX571856.1 L=1550/1585 2217254..2218803 - 16S ribosomal RNA
[12:07:49] Found: 16S_rRNA BX571856.1 L=1551/1585 2039808..2041358 - 16S ribosomal RNA
[12:07:49] Found: 16S_rRNA BX571856.1 L=1548/1585 2335738..2337285 - 16S ribosomal RNA
[12:07:49] Found: 23S_rRNA BX571856.1 L=2919/3232 516231..519149 + 23S ribosomal RNA
[12:07:49] Found: 23S_rRNA BX571856.1 L=2919/3232 561993..564911 + 23S ribosomal RNA
[12:07:49] Found: 23S_rRNA BX571856.1 L=2919/3232 2036412..2039330 - 23S ribosomal RNA
[12:07:49] Found: 23S_rRNA BX571856.1 L=2919/3232 2332376..2335294 - 23S ribosomal RNA
[12:07:49] Found: 23S_rRNA BX571856.1 L=2921/3232 2213965..2216885 - 23S ribosomal RNA
[12:07:49] Found: 5S_rRNA BX571856.1 L=78/119 559007..559084 + 5S ribosomal RNA (partial)
[12:07:49] Found: 5S_rRNA BX571856.1 L=78/119 519226..519303 + 5S ribosomal RNA (partial)
[12:07:49] Found: 5S_rRNA BX571856.1 L=78/119 564988..565065 + 5S ribosomal RNA (partial)
[12:07:49] Found: 5S_rRNA BX571856.1 L=78/119 2036258..2036335 - 5S ribosomal RNA (partial)
[12:07:49] Found: 5S_rRNA BX571856.1 L=78/119 2213810..2213887 - 5S ribosomal RNA (partial)
[12:07:49] Found: 5S_rRNA BX571856.1 L=78/119 2332221..2332298 - 5S ribosomal RNA (partial)
[12:07:49] Found 16 ribosomal RNA features.
[12:07:49] Sorting features and outputting GFF3...
Combined 1 gb file(s)
Results gb file: ./scanned_output/scannedScaffolds.gb
```

Let's try riboSelect:

```
$ riboSelect.py ./scanned_output/scannedScaffolds.gb -o ./select/
2016-12-15 12:15:52 - INFO - Initializing logger
2016-12-15 12:15:52 - INFO - Current usage:
./scanned_output/scannedScaffolds.gb -o ./select/

Loaded 1 records
2016-12-15 12:15:52 - INFO - Processing BX571856.1

2016-12-15 12:15:52 - INFO - hits in BX571856.1

#$ FEATURE rRNA
# Generated cluters for 0 on 20161215
BX571856.1 BX571856.1.fasta1:BX571856.1.fasta2:BX571856.1.fasta3
BX571856.1 BX571856.1.fasta11:BX571856.1.fasta12:BX571856.1.fasta13
BX571856.1 BX571856.1.fasta8:BX571856.1.fasta9:BX571856.1.fasta10
BX571856.1 BX571856.1.fasta14:BX571856.1.fasta15:BX571856.1.fasta16
BX571856.1 BX571856.1.fasta4:BX571856.1.fasta5:BX571856.1.fasta6:BX571856.1.fasta7

```

Now we have our cluster file, reference gb, and data, we can apply `riboSeed.py`…

## riboSeed

```
$ riboSeed.py ./select/riboSelect_grouped_loci.txt -o ./seed/ -F ./SRR4360364_1.fastq -R ./SRR4360364_2.fastq -r ./scanned_output/scannedScaffolds.gb -c 4 -t 1 -v 1
```
