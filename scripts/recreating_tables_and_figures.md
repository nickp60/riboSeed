# Recreating the figures and data used for the manuscript

# Analysis


# Tables


# Figures
##  WGS history plot
Raw data can be downloaded from:
 ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

To generate the plots, navigate to the directory with that file, and run:
`Rscript plot_wgs_history.R prokaryotes.txt`

plot_wgs_history.R can be found under manuscript_results/wgs_history/


## Graphic Abstract
The source files can be found in my lab notebook repo.  If you require access, please contact me (Nick).

## Shannon Entropy within and across genomes
### Sakai
This is pretty straightforward. Download the sakai genome (BA000007.2), and then run riboScan:
```
python3.5 riboScan.py ../../BA000007.2.fasta -o ./sakai_scan/
```
and riboSelect:
```
python3.5  riboSelect.py sakai_scan/scannedScaffolds.gb -o ./sakai_select/
```
and riboSnag:
```
python3.5 riboSnag.py ./sakai_scan/scannedScaffolds.gb ./sakai_select/riboSelect_grouped_loci.txt -o ./sakai_snag_mafft/ --msa_tool mafft -v 1 --title ''
```
### Compare across Sakai genome
The genomes for the comparison were selected using the selectGenomes.R script under manuscript_results/entropy.  In short, this script gets the list of bacterial genomes assemblies from NCBI, selects the non-contig and non-scaffold ones, and gets 25 randomly selected ones (seeded) from the ftp site.

```
Rscript selectGenomes.R assembly_summary.txt ./select_coli_genomes/                                              ```
unzip all the entries in the genomes dir

```
gunzip select_coli_genomes/genomes/
```

Getting the "first" rDNA was a bit more involved:

the gmbH sequence was obtained from Uniprot
```
>gmbH
MAKSVPAIFLDRDGTINVDHGYVHEIDNFEFIDGVIDAMRELKKMGFALVVVTNQSGIA
RGKFTEAQFETLTEWMDWSLADRDVDLDGIYYCPHHPQGSVEEFRQVCDCRKPHPGMLL
SARDYLHIDMAASYMVGDKLEDMQAAVAANVGTKVLVRTGKPITPEAENAADWVLNSLA
DLPQAIKKQQKPAQ
 ```

First, we used a script I wrote to get the coordinates from the gmbH homologue locus, with a similarity threshold for reciprocol blast hits at 95%:
```
python3.5 ~/GitHub/open_utils/orthoML/simpleOrtho.py ./gmbH.faa ./select_coli_genomes/genomes/ -o ortho_gmbH -p 95
```
That outputs a file of the regions in a custom format (ducks under @pjcock's swift left hook) readable by the extractRegion.py utility.


The next stage needs the genomes to be in a single file, so we cat them together:
```
cat ./select_coli_genomes/genomes/* > ./combined_genomes.fa
```

Now, for the tricky part.  The following command reads each set of coordinates from the list outputted by simpleOrtho, finds the sequence id associated with it, and extracts that region with 10kb on either side.  This is because our rDNA of interest is on one side of this gene, but we dont actually know which side, as there may have been inversions, or other funky business.  The extracted regions are written out to a file
```
while read x; do; echo $x | ~/GitHub/open_utils/extractRegion/extractRegion.py ./combined_genomes.fa -e 10000 ; done < ./ortho_gmbH/simpleOrtho_regions.txt  > gmbH_regions.fa
```
Almost there.  So in order to get the rDNA regions back (if they are there), we need to annotate our new set of sequences.  Didn't somebody make a tool to do that?
```
python3.5 riboScan.py ./gmbH_regions.fa -o ./scan_gmbH/ --name gmbH
```
And, group them as operons
```
python3.5 riboSelect.py ./scan_gmbH/scannedScaffolds.gb -o ./select_gmbH/
```
And, extract our rDNAs of interest with a set region flanking them, and make an msa with mafft:
```
python3.5 riboSnag.py ./scan_gmbH/scannedScaffolds.gb ./select_gmbH/riboSelect_grouped_loci.txt -o ./snag_gmbH_mafft/ --msa_tool mafft -v 1 --title ''
```
Thats it!  Our figure is in the output directory.  That wasnt so bad, was it?


### Oh, and the SNP rate for the 25 genomes

First, run Parsnp of the 25 genomes along with the first rDNA from sakai as a reference.
```
ribo scan ./manuscript_results/entropy/sakai_snag_mafft/20170907_region_1_riboSnag.fasta -o manuscript_results/entropy/sakai_1st_rDNA_scan
mkdir manuscript_results/entropy/extraction_gmbH_regions
for i in {1..25} ; do cp manuscript_results/entropy/snag_gmbH_mafft/*_region_${i}_riboSnag.fasta manuscript_results/entropy/extraction_gmbH_regions ; done
 ~/bin/Parsnp-Linux64-v1.2/parsnp -g ./manuscript_results/entropy/sakai_1st_rDNA_scan/scannedScaffolds.gb -d ./manuscript_results/entropy/extraction_gmbH_regions -o 2017-10-13-parsnp
```

Open up the resulting file in Gingr, and export as a vcf

Then, step through getSnpFreqVcf.  The value you are looking for is the "mean(throughout)" one.  Its kind of a hacky way of doing things, but oh well.

## Suppl. Blast results
From the analysis above, where we ran the riboScan, riboSelect, and riboSnag on the E coli sakai genome, the blast results can be used with the script called plotSnagResults.R in the scripts dir.

```
Rscript scripts/plotSnagResults.R ./manuscript_results/entropy/sakai_snag_mafft ./new_blast_results 90
```


## Artificial genome, Representative Mauve output

This is pretty straightforward.


## Running artificial genome analysis
We have written a convenience script called compareColiKleb.sh.  It uses seeded random read generation (with pirs) to make a read set for an  artificial genome created from the rDNAs of E. coli sakai and the surrounding 5kb flaking reigons on either side.  E. coli k12 is then used as a reference to assemble the reads with riboSeed.  We do this 10 times, the results are scored, and then plotted with the script plotSimulatedGenome

```
for i in {1..10}; do ./compareColiKleb.sh $i; done
```


## Running Simulated Genomes for E. coli and Klebsiella

in the scripts folder, you will find two scripts, colisimulation.sh, and klebSimulation.  These orchestrate the running of PIRS and riboSeed.  They are both run as follows:
```
for i in {1..10}; do ./coliSimulation.sh $i; done
```
The results are gathered from the score_reports directory.  we put the combined score reports in folders called "coli" and "kleb", in order to keep track of them.  To generate  the plot, the script plotSimColiKleb.sh was run to extract the relevent data and plot. it is run from RStudio.  Sorry about that...


## Degenerate sequence sequence.

The script that orchestrates this can be found in ./scripts/runDegenerate.sh.  Beware, this is one of the more brain-twisting scripts in the pacakge, so here goes:

It takes a single argument that acts as the replicate number an d also the random seed base (SEED).  The first step is to generate a toy genome with the procedure described previously.  Reads are simulared using PIRS (seeded with SEED).

We empirically determined a good range to test, based off a which gave good resolution for both mutations i the flanking regions and throughout.  This range is (currently) :

FREQS=( 0.0 0.0025 0.0050 0.0075 0.0100 0.0150 0.0200 0.0250 0.0500 0.0750 0.1000 0.1250 0.1500 0.1750 0.2000 0.2250 0.2500 0.2750 0.3000 )

Then, we generate some SUBSEEDS using seedRand.py, one for each of our frequencies.  For each of these, we use riboSim to mutate the flanking regions (or entire sequence) of our toyGenome using riboSim, and concatenating the sequences back together.  This results in a pseudorandomly mutated sequence with a desired frequency.

We then try to run riboSeed to assemble the reads we generated earlier using a progressively worse (more highly mutated)reference sequence.
```
for i in {1..1000}; do ~/GitHub/riboSeed/scripts/runDegenerate.sh $i; done
```
This takes a while (19 frequences, two conditions (ALL, FLANK), 100 replicates, etc).


To generate the plots, use scripts/plotDegenResults:

```
Rscript plotDegenResults.R -r path/to/combinedresults.txt -o ./output_dir/
```

## Degenerate (trusted)
We realized that our intended range overlapped with our cutoff for trusted/untrusted contig treatment of the subseeds.  This required us to run this a second time, but using the old default (--ref_as_contig trust) flag.
