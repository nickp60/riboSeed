# RiboSeed Pipeline
## Description
RiboSeed is an supplemental assembly method to try to address the issue of multiple 16s copies in a genome.  It takes advantage of the fact that while each 16S is identical, the regions flanking are unique, and therfore should be able to be used to seed an alignment.

The pipeline (currently) consists of two main stages:

### Selection and extraction
For generalized selection and extraction of a feature type:

* `otherSelect` [under development]

For ribosomal region selection and extraction

* `riboSelect.py` searches the genome for rRNA annotations, clusters them into likely ribosomal groups, and outputs a colon-separated list of clustered rRNA locus tags by record

NOTE: the format is very simple, and due to the relatively small number of such coding sequences in bacterial genomes, this can be constructed by hand if the clusters do not look appropiate. The format is "genome_sequence_id locus_tag1:locus_tag2", where each line represents a cluster. See example below, where 14 rRNA's are clustered into 6 groups:
```
CM000577.1 FGSG_20052:FGSG_20051:FGSG_20053
CM000577.1 FGSG_20048:FGSG_20047
CM000577.1 FGSG_20049:FGSG_20050
CM000577.1 FGSG_20054:FGSG_20056:FGSG_20055
CM000577.1 FGSG_20058:FGSG_20057
CM000577.1 FGSG_20075:FGSG_20074
```

* `riboSnag.py` takes the list of clustered locus tags and extracts their sequences with flanking regions, optionally turning the coding sequences to N's to minimize bias towards reference. Is used to pull out regions of interest from a Genbank file.  Outputs a directory with a fasta file for each clustered region (and a log file).

* `riboSeed.py` is used to map reads to the extracted regions in an iterative manner, assembling the extracted reads, and then running `SPAdes` assembly to hopefully resolve the contig junctions.

## Known Bugs

* It looks like the `--paired-inference` option will cause an error with `SPAdes` if you are submitting a fastq of singleton/unpaired reads as part of the assembly.

* Submitting `--smalt_scoring` with vastly different scoring schemes usually causes an error.

* `otherSelect` will replace the legacy functioning of `riboSnag` (prior to version 0.8.0)

## Installation
The trickiest part of this whole business is propperly installing SMALT.
Installing with pip3.5 will be the easiest way, but prior to release, clone the repository, and run setup.py.

### Python Requirements:
* Python v3.5 or higher
* Biopython
* pyutilsrnw

### External Requirements
* R (don't ask...)
* SPAdes v3.8 or higher
* SMALT (tested with 0.7.6)
* SAMTools (tested with 1.3.1)
* QUAST (tested with 4.1)

## Note on installation of SMALT
Must have bambamc installed!

#The rest of this readme is old, and will be updated soon.

# riboSnag.py
This script is used to get regions of interest from a Genbank file.  In practice, given a *full* Genbank file (including the sequence), a feature type (rRNA), specific feature type (grepped from '/product='), region (upstream), region length (400bp), this finds all occurrences of a specific feature and extracts the desired sequence, usually a region flanking the feature.

You will probably want to preview your file to figure out the syntax used. (ie, 16s vs 16S, rRNA vs RRNA, etc...)

```
usage: riboSnag.py [-h] [-f FEATURE] [-s SPECIFIC_FEATURES]
                   [-r FEATURE_REGIONS] [-l FEATURE_REGIONS_LENGTHS]
                   [-w WITHIN_FEATURE_LENGTH] [-o OUTPUT]
                   genbank_genome

This is used to extract regions of interest based on coordinates and feature
types.

positional arguments:
  genbank_genome        Genbank file (WITH SEQUENCE)

optional arguments:
  -h, --help            show this help message and exit
  -f FEATURE, --feature FEATURE
                        Feature, such as CDS,tRNA, rRNA
  -s SPECIFIC_FEATURES, --specific_features SPECIFIC_FEATURES
                        colon:separated -- specific features to be grepped
                        from product, such as 16S or tRNA-Ala
  -r FEATURE_REGIONS, --feature_regions FEATURE_REGIONS
                        colon:separated -- upstream, upstream_within,
                        downstream, downstream_within, identity, identityN
  -l FEATURE_REGIONS_LENGTHS, --feature_regions_lengths FEATURE_REGIONS_LENGTHS
                        colon:separated -- length of region
  -w WITHIN_FEATURE_LENGTH, --within_feature_length WITHIN_FEATURE_LENGTH
                        if using up/downstream_within, where you return some
                        of the feature itself, this dictates how many bp
  -o OUTPUT, --output OUTPUT
                        output directory

```

## Output
Output is a directory with a fasta file for each of the regions found in the genome with a header stating the coordinates of the feature and the coordinates of the sequence extracted.


# riboSeed.py
This takes the output from `riboSnag`, and for each file it maps reads to it, extracts mapped reads, assembles extracted reads, and then repeats that depending on the number of iterations specified.  Then, it uses `SPAdes` to assemble the raw reads with the other fragments as "trusted contigs".

```

usage: riboSeed.py [-h] [-F FASTQ1] [-R FASTQ2] [-S FASTQS] [-n EXP_NAME]
                   [-m METHOD] [-c CORES] [-r REFERENCE_GENOME] [-o OUTPUT]
                   [--paired_inference] [--subtract] [--keep_unmapped]
                   [--ref_as_contig REF_AS_CONTIG] [--temps] [-i ITERATIONS]
                   [-v VERBOSITY] [--DEBUG] [--force] [--noclobber]
                   [--smalt_scoring SMALT_SCORING] [--spades_exe SPADES_EXE]
                   [--samtools_exe SAMTOOLS_EXE] [--smalt_exe SMALT_EXE]
                   [--quast_exe QUAST_EXE]
                   seed_dir

Given regions from riboSnag, assembles the mapped reads

positional arguments:
  seed_dir              path to roboSnag results directory

optional arguments:
  -h, --help            show this help message and exit
  -F FASTQ1, --fastq1 FASTQ1
                        forward fastq reads, can be compressed
  -R FASTQ2, --fastq2 FASTQ2
                        reverse fastq reads, can be compressed
  -S FASTQS, --fastq_single FASTQS
                        single fastq reads
  -n EXP_NAME, --experiment_name EXP_NAME
                        prefix for results files; default: riboSeed
  -m METHOD, --method_for_map METHOD
                        availible mappers: smalt; default: smalt
  -c CORES, --cores CORES
                        cores for multiprocessing workers; default: 1
  -r REFERENCE_GENOME, --reference_genome REFERENCE_GENOME
                        fasta reference genome, used for estimating insert
                        sizes, QUAST, and SPAdes
  -o OUTPUT, --output OUTPUT
                        output directory; default: /home/nicholas/GitHub/FB/Ec
                        oli_comparative_genomics/scripts
  --paired_inference    if --paired_inference, mapped read's pairs are
                        included; default: False
  --subtract            if --subtract, reads aligned to each reference will
                        not be aligned to future iterations. Probably you
                        shouldnt do thisunless you really happen to want to
  --keep_unmapped       if --keep_unmapped fastqs are generated containing the
                        unmapped reads; default: False
  --ref_as_contig REF_AS_CONTIG
                        if 'trusted', SPAdes will use the seed sequences as a
                        --trusted-contig; if 'untrusted', SPAdes will treat as
                        --untrusted-contig. if '', seeds will not be used
                        during assembly. See SPAdes docs; default:
  --temps               if --temps, intermediate files will be kept; default:
                        False
  -i ITERATIONS, --iterations ITERATIONS
                        if iterations>1, multiple seedings will occur after
                        assembly of seed regions ; default: 2
  -v VERBOSITY, --verbosity VERBOSITY
                        1 = debug(), 2 = info(), 3 = warning(), 4 = error()
                        and 5 = critical(); default: 2
  --DEBUG               if --DEBUG, test data will be used; default: False
  --force               if --force, existing results dirs will be used;
                        default: False
  --noclobber           if --noclobber, results dirs will be overwritten, not
                        deleted and written fresh; default: False
  --smalt_scoring SMALT_SCORING
                        submit custom smalt scoring via the smalt -S scorespec
                        option; default: match=1,subst=-4,gapopen=-4,gapext=-3
  --spades_exe SPADES_EXE
                        Path to spades executable; default: spades.py
  --samtools_exe SAMTOOLS_EXE
                        Path to bwa executable; default: samtools
  --smalt_exe SMALT_EXE
                        Path to smalt executable; default: smalt
  --quast_exe QUAST_EXE
                        Path to quast executable; default: quast.py
```


## Output

This outputs two main directories: `map` and `results`.  If `--temps` is true, temporary files from the mapping scheme will be retained, and is useful for assessing problems.
The results directory will contain a 'mauve' directory with all the extended fragments, the mapped fastq files, and a `de_novo` and `de_fere_novo` folder, containing the results with the *de novo* mapping and supplemented mapping, respectively.


## Suggested Running
There are a lot of commandline options for these, so it can help to run concurrently as  script.

```

ITER=2
LENGTH=1000
WIDTH=1000
SEED_DIR="/home/nwaters/coli_w_1000_l_1000
OUT_DIR="/home/nwaters/Seed_w_1000_l_1000_i_2

python3.5 riboSnag.py /home/nwaters/GitHub/FB/Ecoli_comparative_genomics/scripts/riboSeed_pipeline/NC_011751.1.gb -f rRNA -s 16S:5S -r upstream_within:downstream_within -l ${LENGTH}:${LENGTH} -w ${WIDTH} -o ${SEED_DIR}

python3.5 riboSeed.py -F ~/GitHub/FB/Ecoli_comparative_genomics/scripts/riboSeed_pipeline/WTCHG_142097_201102_1.trimmed.fq.gz -R ~/GitHub/FB/Ecoli_comparative_genomics/scripts/riboSeed_pipeline/WTCHG_142097_201102_2.trimmed.fq.gz -S ~/GitHub/FB/Ecoli_comparative_genomics/scripts/riboSeed_pipeline/WTCHG_142097_201102_singles.trimmed.fq.gz -n two_iter -m smalt -c 12 -o ${OUT_DIR} ${SEED_DIR}  --reference_genome ~/GitHub/FB/Ecoli_comparative_genomics/scripts/riboSeed_pipeline/NC_011751.1.fasta --ref_as_contig trusted  --iterations ${ITER} -v 1 --smalt_scoring match=1,subst=-5,gapopen=-5,gapext=-3

```
