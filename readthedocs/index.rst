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

Preprint of the riboSeed manuscript can be found
`here <http://www.biorxiv.org/content/early/2017/07/14/159798>`__;
comments welcome!

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
