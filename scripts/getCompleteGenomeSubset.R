#  get genomes
version = 0.2

# From the list of sequenced genomes from the NCBI FTP site, extract ones to search for the RBHs from
#install.packages("data.table")
library(data.table)

print("USAGE: Rscript selectGenomes.R assembly_summary.txt path.to/output/dir/ genus species")
nstrains <- 25
args <- commandArgs(trailingOnly = TRUE)
print(args)
assembly_summary <- as.character(args[1])
outdir <- file.path(getwd(), as.character(args[2]))
genus <- as.character(args[3])
species <- as.character(args[4])
dir.create(outdir, showWarnings = T)
dir.create(file.path(outdir, "genomes"), showWarnings = T)


search_term<-"Escherichia coli"
search_term  <- paste(genus, species)
get_from_R <- T
# setwd("~/")
# outdir<-"~/GitHub/riboSeed/manuscript_results/entropy/coli/"
# system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt")
raw_summaries<-fread(assembly_summary, skip = 1)


full<-raw_summaries[raw_summaries$genome_rep == "Full"]

complete<- full[!full$assembly_level %in% c("Contig", "Scaffold"),]

not_excluded<-complete[complete$excluded_from_refseq == "",] 

of_interest <-not_excluded[grepl(search_term, not_excluded$organism_name),]

write.table(paste0(of_interest$ftp_path, "_genomic.fna.gz"),
            file = file.path(outdir, paste0(gsub(" ","_", search_term),"_genomes_of_interest.txt")), 
            row.names = F, col.names = F, quote = F)

set.seed(27)
if (nrow(of_interest) > nstrains){
  subset_of_interest <- of_interest[sample(nrow(of_interest), nstrains), ]
} else {
  subset_of_interest <- of_interest
}

if (get_from_R){ 
  for (i in subset_of_interest$ftp_path){
    file <- paste0(i, "/", gsub("(.*)/(.*$)", "\\2", i), "_genomic.fna.gz") # or _genomic.gbff.gz 
    outfile <- file.path(outdir, "genomes", paste0(gsub("(.*)/(.*$)", "\\2", i), "_genomic.fna.gz")) # or _genomic.gbff.gz 
    system(paste0("wget ",file, " -O ", outfile))
  }
} else {
  subset_of_interest$ftp_path_full<- paste0(of_interest$ftp_path, "/", gsub("(.*)/(.*$)", "\\2", of_interest$ftp_path), "_genomic.fna.gz") # or _genomic.gbff.gz 
  write.table(subset_of_interest$ftp_path_full,
              file = paste0(outdir, "/", gsub(" ","_", search_term),"_genomes_of_interest_fullpath.txt"), 
              row.names = F, col.names = F, quote = F)
  
}

