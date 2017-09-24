# This script is run to recreate the subset of sthe Prokaryotes.txt file used
# by cgfind.  

proks <- read.csv("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt", header=T, sep="\t", stringsAsFactors=FALSE)
proks <- proks[grepl("Complete", proks$Status), ]

# proks_repr <- read.csv2("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prok_representative_genomes.txt", header=T, sep="\t", stringsAsFactors=FALSE)
# proks_repr <- proks_repr[proks_repr$X.Species.genus != "", ]  

proks$genus <- gsub("(.*?) (.*?) (.*)", "\\2", proks$X.Organism.Name)
proks$strain <- gsub("(.*?) (.*?) (.*)", "\\1", proks$X.Organism.Name)
proks$acc <- gsub("(.*?)/.(.*)", "\\1", 
     gsub("chromosome:", "", proks$Replicons))

proks$full_ftp <- paste0(proks$FTP.Path, "/",
                         basename(proks$FTP.Path), "_genomic.fasta")


for (thing in c("Status", "WGS", "Group", "SubGroup", "BioProject.Accession", "BioSample.Accession", "BioProject.ID", "Reference", "Pubmed.ID", "Center", "Proteins", "Strain", "TaxID", "FTP.Path", "Replicons")){
  proks[ , thing] <- NULL
  print(object.size(proks))
}
write.csv(x = proks, file = paste0(Sys.Date(), "_complete_genomes_min.csv"))
