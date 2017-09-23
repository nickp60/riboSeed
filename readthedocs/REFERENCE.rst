Choosing an Appropriate Reference
===============

Introduction
----

Any sort of analysis that involves a reference strain begs the question: which strain should I use? Luckily, riboSeed's *de fere novo* approach gives you a bit of flexiblility here, in that the method is not affected my major genomic restructuring events, as long as they do not occur in the rDNA region.  Still, using a close reference maximizes your chance for a successful assembly.

There are two ways that we recommend selecting a reference.  One is a more robust method that provides a high degree of certainty, but requires a bit of legwork. The other method can be done entirely through your web browser, but is much less robust. We recommend method #1, but if you have already pulled out all of your hair from bioinformatics software installation, method #2 is a nice, painless alternative.

Method #1: Kraken
~~~


Kraken is a kmer-based phylogeny tool that can be used to idenify the strains present in a metagenomic dataset;  the installation and usage `instructions can be found here __ <https://ccb.jhu.edu/software/kraken/>`__

- Download and install Kraken, along with the MiniKraken database from their website.
- Run Kraken on your isolate's reads, and generate the Kraken report.

::

    kraken --db MiniKraken reads1.fq reads2.fastq  > sequences.kraken
    kraken-translate --db MiniKraken sequences.kraken > sequences.labels
    kraken-report --db MiniKraken sequence.kraken

Because the MiniKraken database was built from all the complete genomes from RefSeq, it should be easy to identify which strain in the database has the closest match to your sequenced isolate.

PS:  This is a great time to check if you have any contamination in your sample;  thanks, Kraken!


Method #2: reads2type and cgfind
~~~


reads2type is also a kmer-based phylogeny tool, but it relies on a lightweight, prebuilt database that allows the analysis to be performed in your web browser, and it doesn't require you to upload your whole read file to a webserver.  It works by taking one read at a time from your file, generating 55-mers, and comparing to its databse. If there is not enough taxanomic information to indentify the isolate off of that read alone, additional reads will be processed until a single taxonomy is achieved.  This method works best on trimmed reads. `Instructions and the webserver can be found here __ <https://cge.cbs.dtu.dk/services/Reads2Type/>`__

Once you have a genus and species, you can use ``cgfind``, a tool we developed to provide easy access to downloadable genomes based on the complete prokaryotic genomes found in NCBI.  `it can be found here __<https://cge.cbs.dtu.dk/services/Reads2Type/>`__  Just enter your genus and species name, and select one of the avaible strains to download.
