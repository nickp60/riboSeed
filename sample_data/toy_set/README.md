# Toy dataset example
This is a example dataset used to show the utility of the riboSeed pipeline.  We extracted all the rRNA operons from the E. Coli Sakai genome (BA000007.2) with 5000bp flanking regions.  Then, we manually combined those regions together into a single fasta file, our pseudogenome.

Then, we started the pipeline.  The `scanScaffolds.sh` preprocessing script uses Barrnap to annotate our psuedogenome with rRNA's, and used `riboSelect.py` to generate the cluster file
```~/GitHub/riboSeed/scanScaffolds.sh ./pseudata/ fa ./pseudata/scanned/ bac ```
```python3.5 ~/GitHub/riboSeed/riboSeed/riboSnag.py ./pseudata/scanned/scanScaffolds_combined.gb ./pseudata/select/20161005_riboSelect_grouped_loci.txt --clobber -v 1 -o ./pseudata/snag/```

We extracted those regions and the default 700bp flanking regions with riboSnag.
```python3.5 ~/GitHub/riboSeed/riboSeed/riboSelect.py ./pseudata/scanned/scanScaffolds_combined.gb -o ./pseudata/select/ --clobber -v 1```
We used ART to generate synthetic reads (145bp paired-end reads with an average 500bp insert size, to a coverage depth of 50x.
```~/bin/art_bin_MountRainier/art_illumina -ss MSv3 -i ./pseudata/combined_colapsed.fa -l 145 -p -m 500 -s 20 -f 50 -o toy_reads```
Riboseed was run on the dataset, and then visualized with mauve.
```python3.5 ~/GitHub/riboSeed/riboSeed/riboSeed.py  ./pseudata/snag/ -F ./toy_reads1.fq -R toy_reads2.fq -r ./data/littlegenome.fa -v 1 -i 3 --ref_as_contig trusted -c 4 -o ./pseudata/seed1/```


```
~/GitHub/riboSeed/scanScaffolds.sh ./data/seq/ .fasta ./data/scanned/ bac
~/GitHub/riboSeed/riboSeed/riboSelect.py ./data/scanned/scanScaffolds_combined.gb -o ./data/select/
python3.5 ~/GitHub/riboSeed/riboSeed/riboSnag.py ./data/scanned/scanScaffolds_combined.gb ./data/select/20161005_riboSelect_grouped_loci.txt -l 5000 -o ./data/seed/
cat ./data/seed/20161005_region_* >./data/littlegenome.fa
\# manually smooshed the contigs together by deleting headers and up to 50 bp sequence to preserve line lengths, called it combined_collapsed, and put it in ./pseudata/



├── data  # Starting genome from E. coli Sakai
│   ├── BA000007.2.gb
│   ├── combined_colapsed.fa
│   ├── combined.fa
│   ├── littlegenome.fa
│   ├── scanned
│   ├── seed
│   ├── select
│   ├── seq
│   └── snag
├── pseudata  # using the combined_collapsed.fa genone
│   ├── combined_colapsed.fa
│   ├── scanned
│   ├── seed1
│   ├── select
│   └── snag
├── README # this file
├── toy_reads1.aln
├── toy_reads1.fq
├── toy_reads2.aln
└── toy_reads2.fq
```
