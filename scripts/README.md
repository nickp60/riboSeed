# Recreating the figures and data used for the manuscript

<!-- # Analysis -->


<!-- # Tables -->


<!-- # Figures -->
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
The genomes for the comparison were selected using the selectGenomes.R script (renamed getCompleteGenomeSubset.R, under "scripts" dir) under manuscript_results/entropy.  In short, this script gets the list of bacterial genomes assemblies from NCBI, selects the non-contig and non-scaffold ones, and gets 25 randomly selected ones (seeded) from the ftp site.

```
Rscript selectGenomes.R assembly_summary.txt ./select_coli_genomes/
```

unzip all the entries in the genomes dir

```
gunzip select_coli_genomes/genomes/
```

Getting the "first" rDNA was a bit more involved:

Manual inspection of the genome showed the first rDNA had a neighboring gene called gmbH. A gmbH sequence was obtained from Uniprot
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
We realized that our intended range overlapped with our cutoff for trusted/untrusted contig treatment of the subseeds.  This required us to run this a second time, but using the old default (--ref_as_contig trusted) flag.


## Pseudomonas
### Assembly
riboSeed was run on the dataset as follows:
```
/ribo run -r ./data/NZ_CP017149.1.fasta -o /home/nick/results/2017-10-25-pao_ATCC_ref/ -F ./data/SRR3500543_1.fastq -R ./data/SRR3500543_2.fastq --memory 42 --cores 12 --threads 2
```

To determine how it performed in relation to using the whole genome as a reference, we reran the last SPAdes call as follows:

```
/python3.5 /home/nick/miniconda3/bin/spades.py -t 6 -m 21 --careful -k 21,33,55,77 --pe1-1 ./data/SRR3500543_1.fastq --pe1-2 ./data/SRR3500543_2.fastq --trusted-contigs /home/nick/results/2017-10-25-pao_ATCC_ref/scan/contigs/NZ_CP017149.1_0.fasta  -o /home/nick/results/2017-10-25-pao_ATCC_ref/seed/final_full_ref_assembly
```
### QUAST
Then, after copying the reesults to the mauve dir, we ran quast on the results as follows:

```
/home/nicholas/bin/quast-4.4/quast.py -o 2017-10-25-pao_ATCC_ref/QUAST/  -R ./CP015377.1.fasta --threads 4 ./2017-10-25-pao_ATCC_ref/NZ_CP017149.1_de_fere_novo_contigs.fasta ./2017-10-25-pao_ATCC_ref/NZ_CP017149.1_de_novo_contigs.fasta ./2017-10-25-pao_ATCC_ref/NZ_CP017149.1_full_ref.fasta
```

The table was then recapitualed in LaTeX.  Bummer, right?  I wish they could do colored output in latex by default.


### SNPs
To count the snps comparing the reference to the de fere novo assembly and the de novo assembly, we ran parsnp as follows:

```
~/bin/Parsnp-Linux64-v1.2/parsnp -g ./CP015377.gb -c -d ./2017-10-25-pao_ATCC_ref/parsnp_genomes/ -p 8 -o ./2017-10-25-pao_ATCC_ref/parsnp_results/
```
where the parsnp_genomes dir contained the de fere novo assembly and the refernece as a fasta.

The resulting file was opened in gingr and exported as a vcf.

riboScan was run on the true reference (CP015377.1.fasta) to generate a gff file, which is then used along with the vcf in a scrript poorly named checkParsnpVcf.R.

## Verifying B. cereus vd118
### The problem.
We noticed that riboSeed calculated the percent mapped reads at less than 80% for this strain.

### Kraken
Upon further investigation with Kraken, we determined that 30% of the reads were unclassifiable with the Kraken miniKraken database, and only 70% of the trimmed reads were classifiable to bacteria.  Of those, 68% appeared to be B. cereus.  Thats not quite enough for my taste.

```
kraken --db ~/bin/minikraken_20141208/ ~/Downloads/hi/cereus/trimmed/insert_180_1__cov250x.fastq ~/Downloads/hi/cereus/trimmed/insert_180_1__cov250x.fastq > ~/Downloads/hi/cereus/2017-10-18-minikraken_cereus  _README
kraken --db ~/bin/minikraken_20141208/ ~/Downloads/hi/cereus/raw/insert_180_1__cov250x.fastq ~/Downloads/hi/cereus/raw/insert_180_2__cov250x.fastq > ~/Downloads/hi/cereus/2017-10-18-minikraken_cereus_untrimmed  _README
kraken-report --db ~/bin/minikraken/ ~/Downloads/hi/cereus/2017-10-18-minikraken_cereus > ~/Downloads/hi/cereus/2017-10-18-minikraken_cereus_report
kraken-report --db ~/bin/minikraken/ ~/Downloads/hi/cereus/2017-10-18-minikraken_cereus_untrimmed > ~/Downloads/hi/cereus/2017-10-18-minikraken_cereus_untrimmed_report
```

### metaspades
We assembled the genome denovo with metaspades for downstream analysis with MBBC, MaxBin, and blobtools

```
metaspades.py -o ./2017-10-17-metaspades-cereus --pe1-1 ~/Downloads/hi/cereus/trimmed/insert_180_1__cov250x.fastq --pe1-2 ~/Downloads/hi/cereus/trimmed/insert_180_2__cov250x.fastq --threads 12 --memory 50
```

### MAxBin

### blobtools
We wanted to generate a blobplot with blobtools.  That requires 1) an assembly 2) a mapping file, and 3) a "hits file.  After much jockeying, I eventualy got thisworking after downloading the entire nt database to microgue.

We genereated the blast hits file with the dollowing command:

```
blastn -db ~/BLASTDB/nt -query ~/results/2017-10-17-metaspades-cereus/contigs.fasta -out results/2017-10-19-cereus_assembly_blastn_nt_ids.out -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -num_threads 12
```
And then we rna blobtools as follows
```
# build the database
blobtools create -i ~/results/2017-10-17-metaspades-cereus/contigs.fasta -y spades -t ~/results/2017-10-19-cereus_assembly_blastn_nt_ids.out -o ./blobplot_spades_cov_tax
# generate the table (for some reason? just for viewing, I guess?)
grep '^##' blobplot_spades_cov_tax.blobDB.table.txt > blobplot_spades_cov_tax_view ;  grep -v '^##' blobplot_spades_cov_tax.blobDB.table.txt |  column -t -s $'\t'  >> blobplot_spades_cov_tax_view
# build the blobplot at species level
blobtools blobplot  -i blobplot_spades_cov_tax.blobDB.json -r species  -o ./
```





### Maxbin
<!-- Now that we had proof of contamination, we used MAxbin to separated the contigs likely belonging to each assembly -->


We attempted to use maxbn to separate the reads belonging to each of the groups, but presumably because the GC content was so similar and the fact that even after filtering out the reads directly mapping to the B. cereus reference, many of the resulting contigs showed homology to B. cereus, we decided this was not the most appropriate way to isolate the contmainates.

### Mapping to filter out Non-B. cereus reads:
We were able to make it a bit clearer by first mapping the reads to the reference genome, filtering out those reads mapping, and assembling the remaining reads.
```
bwa index ~/Downloads/hi/cereus/AE017194.1.fasta
bwa mem ~/Downloads/hi/cereus/AE017194.1.fasta ~/Downloads/hi/cereus/trimmed/insert_180_1__cov250x.fastq ~/Downloads/hi/cereus/trimmed/insert_180_2__cov250x.fastq > cereus_subset_atcc.sam
samtools fastq -F 12 cereus_subset_atcc.sam -1 cereus_subset_atcc_mappedF12_1.fastq -2 cereus_subset_atcc_mappedF12_2.fastq

metaspades.py --pe1-1 ./cereus_subset_atcc_mappedF12_1.fastq --pe1-2 ./cereus_subset_atcc_mappedF12_2.fastq -o ./assembled_atcc/
```
Then, we reran the blast search toget a hits file for blobtools, and generated a blobplot

```
blastn -db ~/BLASTDB/nt -query ./assembled_atcc/contigs.fasta -out atcc_blastn_nt_ids.out -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -num_threads 12
blobtools create -i ./assembled_atcc/contigs.fasta  -y spades -o ./blob_atcc/ -t ./atcc_blastn_nt_ids.out
```


That was all a bit confusing, so we made a script (./scripts/blob_cereus.sh)that runs all the commands.  Paths are hardcoded, so it is more of a description of the work than a portable script, but oh well.


## Archaea
### Formicum
This was a strange one:  we found out after enforcing stricter kmer choices that performace dropped off.
riboSeed was run, and then we re-ran spades with different kmers
```
 ribo seed -r ./scan/scannedScaffolds.gb -S1 ./DRR017790_1.fastq  -z --cores 4 --memory 16 -o form_infer --ref_as_contig infer ./select/riboSelect_grouped_loci.txt
/home/nw42839/.pyenv/shims/python3.5 /home/nw42839/miniconda3/bin/spades.py -t 16 -m 32 --careful -k 21,33,55,77,99 --pe1-s ./DRR017790_1.fastq --trusted-contigs /home/nw42839/Results/form/form_infer_k85/final_long_reads/riboSeedContigs.fasta  -o /home/nw42839/Results/form/form_infer_k85/final_de_fere_k99

```



### Barkeri
This dataset is HUGE, so thats why we downsampled it to just 5\% of the reads:
```
~/miniconda3/bin/fastq-dump --split-files SRR2064286

seqtk sample -s 27 ./SRR2064286_1.fastq 0.05 > SRR2064286_1_sub.fastq; seqtk sample -s 27 ./SRR2064286_2.fastq 0.05 > SRR2064286_2_sub.fastq
```
Then, we ran riboSeed:

```
/home/nw42839/miniconda3/bin/ribo scan CP009526.1.fasta -o scan
/home/nw42839/miniconda3/bin/ribo select ./scan/scannedScaffolds.gb -o select
/home/nw42839/miniconda3/bin/ribo seed ./select/riboSelect_grouped_loci.txt -r ./scan/scannedScaffolds.gb -F ./SRR2064286_1_sub.fastq -R ./SRR2064286_2_sub.fastq --ref_as_contig trusted --cores 4 --memory 24 -z -o seed_trusted
```


## More entropy the suppl figures:
for generating the entropy figures within and across genomes, we updated things to use a script: runEntropyComparisons.sh.  Its a clunky one, but it works.  It reads the entropy_manifest.tab file, and requires a combined archaeal/bacterial assemblies_summary.txt file from NCBI.  The beginning of the script has comments about how to regenerate that file.  In short, I manually inspected each genome and picked a gene close the the first rDNA operon.  The script picks up to 25 randome genomes of the same genus and species, and runs the analysis outlined above (extract those regions + 10bk up and down, find the rDNAs, extract the rDNAs, allign with mafft, plot).  It also plots the within-genome rDNA entropy for that genome.  We got the sense that the supplementary material was too short, and figured this would lengthen things substantially.  (Just kidding, a reviewer asked).

## The ribo structure analysis
We found several genomes with strange  rDNA structures, and downloaded them manually from ncbi.  Then, after putting them in a folder called "./2018-01-26-odd_rDNAs/", we ran riboStrcuture:

```
ribo structure ./2018-01-26-odd_rDNAs/ -o oddballs
```

as a comparison, because we had the genomes already for the entropy experiment, we copied those all to a new dir:

```
cp 2018-01-17-within_vs_across/*/ref/*.fasta ./2018-01-27-normal/
```

and ran ribostructure


```
rm normal/ -r; ribo structure ./2018-01-27-normal/ -o normal
```
