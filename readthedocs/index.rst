.. riboSeed documentation master file, created by
   sphinx-quickstart on Wed Aug 16 23:19:52 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to riboSeed's documentation!
====================================
                 |ribologo|

Impatient? See our `Quickstart Guide <https://github.com/nickp60/riboSeed/blob/master/quickstart.md>`__

A brief overview of the theory can be found
`here <https://nickp60.github.io/riboSeed.html>`__.

The manuscript can be found here:
`Nicholas R Waters, Florence Abram, Fiona Brennan, Ashleigh Holmes, Leighton Pritchard; riboSeed: leveraging prokaryotic genomic architecture to assemble across ribosomal regions, Nucleic Acids Research, Volume 46, Issue 11, 20 June 2018, Pages e68, https://doi.org/10.1093/nar/gky212 <https://academic.oup.com/nar/article/46/11/e68/4955760>`__

Description
-----------

riboSeed is an supplemental assembly refinement method to try to address
the issue of multiple ribosomal regions in a genome, as these create
repeates unresolvable by short read sequencing. It takes advantage of
the fact that while each region is identical, the regions flanking are
unique, and therefore can potentially be used to seed an assembly in
such a way that rDNA regions are bridged.


Contents:

.. toctree::
   :maxdepth: 4

   OVERVIEW
   INSTALLATION
   REFERENCE
   USAGE
   PIPELINE
   ACCESSORY


.. automodule::
   :members:
   :show-inheritance:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |ribologo| image:: logo_1.svg
   :target: https://github.com/nickp60/riboSeed/blob/master/icon/logo_1.svg
