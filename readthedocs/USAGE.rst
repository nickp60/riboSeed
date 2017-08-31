Usage
===============

Minimal example: ``run_riboSeed.py``
-------

The pseudogenome was constructed from the 7 rDNAs separated by several kb of flanking DNA.  If can be found under `./riboSeed/integration_test/concatenated_seq.fasta`.  If you have installed using setuptools, the `integration_test` folder will be installed in the site-packages dir, such as `/venv-riboSeed/lib/python3.5/site-packages/riboSeed/integration_data/`.

Two read files can be found in the same directory.

To run the whole riboSeed pipline, use the following command:

::

    python3.5 ./riboSeed/run_riboSeed.py ./riboSeed/integration_data/concatenated_seq.fasta -F ./riboSeed/integration_data/test_reads1.fq -R ./riboSeed/integration_data/test_reads2.fq -o ./test1/ -v 1


Whats going on:
~~~~
``run_riboSeed.py`` is used to run the pipeline with the most commonly used settings. It first creates a config file, tracking down your system executables
for the required tools, and setting the default parameters for things not
specified as args to run_riboSeed.

Then, riboScan is run to re-annotate your reference, riboSelect calls the rDNA
operons, and riboSeed runs the *de fere novo* assembly.

If you want to change the behaviour of the programs under the hood, all of the
command line options not set by run_riboSeed are defined in the config file in
the output directory. After editing the parameters in the config file, you can
submit it to run_riboSeed using the -c flag.

Running individual scripts
-------

All of the scripts in the package can be run individually. Perhaps you want to
modify barrnap's behaviour in riboScan, or you want to experiment with
different feature selectors in riboSelect.  Go for it!
