Description
===============

riboSeed is an supplemental assembly refinement method to try to address
the issue of multiple ribosomal regions in a genome, as these create
repeates unresolvable by short read sequencing. It takes advantage of
the fact that while each region is identical, the regions flanking are
unique, and therefore can potentially be used to seed an assembly in
such a way that rDNA regions are bridged.

The Problem
---------------
As you probably know, repeated regions are difficult to resolve when sequencing with a short read technology. Specifically, if the length of the repeat exceeds the length of the kmers used to construct the de Bruijn graph, the repeat cannot be resolved.

rDNAs, the genomic regions containing the sequences coding for ribosomal RNAs, are often found multiple times in a single genome. The rDNAs are usually around 5kb long, which is much longer than the length of short reads. Because of this, these regions cause breaks in the assembly.

Due to how well rDNA is conserved within a taxa, we hypothesized that if the regions flanking the rDNAs are sufficiently unique within a genome, those regions would be able to locate an rDNA within the genome during assembly.

The Method
---------------

We call our method a *de fere novo* assembly (meaning "starting from almost nothing"), as we use a subassembly technique to minimize the bias caused by reference choice.  We map the short reads to the reference genome, extract the reads mapping to rDNA (with flanking) regions, and perform subassemblies with SPAdes to reassemble the rDNA and flanking regions from the reads.  These "long reads" are concatenated together separated with 5kb of N's. The reads are then mapped to the concatenated sequence and and subassembled for several additional iterations.


The Results
---------------

We generated a simulated genome from the 7 rDNA regions with 5kb flanking regions, and then used [ART (MountRainier-2016-06-05)]("https://www.niehs.nih.gov/research/resources/software/biostatistics/art/") to generated simulated MiSeq reads of various depths.

|sim|

In this ```Mauve``<http://darlinglab.org/mauve/mauve.html>`__ visualization, we show (from top to bottom) the reference simulated genome, riboSeed's  *de fere novo* assembly,  *de novo* assembly, and a negative control  *de fere novo* assembly using a *Klebsiella* reference genome.  The results show that with riboSeed's *de fere novo* assembly correctly joins six of the seven rDNA regions to reconstruct the simulated genome with only short reads.  By contrast, the short reads alone failed to bridge any gaps caused by the repeated rDNAs, and the assembly using a poor reference choice only assembled across a single rDNA region. We have run this successfully on many real datasets with positive results


.. |sim| image:: mauve_simulated.png
   :target: https://github.com/nickp60/riboSeed/blob/master/docs/images/mauve_simulated.png
