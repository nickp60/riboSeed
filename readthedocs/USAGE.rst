Usage
===============

Minimal example: ``ribo run``
-------

The pseudogenome was constructed from the 7 rDNAs separated by several kb of flanking DNA.  If can be found under `./riboSeed/integration_test/concatenated_seq.fasta`.  If you have installed using setuptools, the `integration_test` folder will be installed in the site-packages dir, such as `/venv-riboSeed/lib/python3.5/site-packages/riboSeed/integration_data/`.

Two read files can be found in the same directory.

To run the whole riboSeed pipline, use the following command:

::

    ribo run ./riboSeed/integration_data/concatenated_seq.fasta \
              -F ./riboSeed/integration_data/test_reads1.fq \
              -R ./riboSeed/integration_data/test_reads2.fq \
              -o ./test1/ -v 1


Whats going on:
~~~~
``ribo run`` is used to run the pipeline with the most commonly used settings. It first creates a config file, tracking down your system executables
for the required tools, and setting the default parameters for things not
specified as args to run_riboSeed.

Then, ``ribo scan`` is run to re-annotate your reference, ``ribo select`` calls the rDNA
operons, and ``ribo seed`` runs the *de fere novo* assembly.

If you want to change the behaviour of the programs under the hood, all of the
command line options not set by ``ribo run`` are defined in the config file in
the output directory. After editing the parameters in the config file, you can
submit it to ``ribo run`` using the -c flag.

Running individual scripts
-------

All of the elements of the package can be run individually: Perhaps you want to
modify barrnap's behaviour in ``scan``, or you want to experiment with
different feature selectors in ``select``.  Go for it!

::

   $ ribo

   Description: A suite of tools to perform de fere novo assembly to bridge
   gaps caused by rDNA repeats

   Usage:   ribo <command> [options]

   Available commands:
   -run        execute pipeline (scan, select, seed, sketch, and score)
   -scan       reannotate rRNAs in a FASTA file
   -select     group rRNA annotations into rDNA operons
   -seed       perform de fere novo assembly
   -snag       extract rDNA regions and plot entropy
   -sim        perform simulations used in manuscript
   -sketch     plot results from a de fere novo assembly
   -stack      compare coverage depth in rDNA regions to rest of genome
   -score      score batches of assemblies with BLASTn
   -swap       swap contigs from assemblies
   -spec       use assembly graph to speculate number of rDNAs
   -structure  view the rRNA operon structure across several genomes
   -config     write out a blank config file to be used with `run`
