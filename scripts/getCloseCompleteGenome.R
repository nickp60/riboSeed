# This script is used to get the closes complete genome based on some very basic heuristics

help <- ' A basic taxid to lineage flatfile can be found here
https://gitlab.com/zyxue/ncbitax2lin-lineages/blob/master/lineages-2017-03-17.csv.gz

#Download and unpack it, read it in
db <- read.csv2("~/Downloads/lineages-2017-03-17.csv", header = T, stringsAsFactors = F, sep=",")

# read in the prokaryotes.csv from NCBI
raw <- read.csv2("~/GitHub/riboSeed/manuscript_results/wgs_history/prokaryotes.txt", header = T, stringsAsFactors = F, sep="\t")

# tidy up
raw <- raw[raw$Release.Date != "-",]
raw <- raw[grepl("Complete", raw$Status), ]

# merge, and subset columns

annotated <- merge(db, raw, by.x = "tax_id", by.y = "TaxID", all.x = F, all.y = T)  


annotated <- annotated[!is.na(annotated$superkingdom), 
                       c(
                         colnames(annotated)[1:8], # taxid and taxonmy superkingdom through species 
                        "BioProject.Accession",
                         "Size..Mb.",
                         "GC.", "Replicons", "Scaffolds", 
                         "Genes", "Proteins", "Release.Date", "Modify.Date",
                        "Status", "Center", "BioSample.Accession", "Assembly.Accession", "Reference", "FTP.Path")]

# and finally, write out the combined all_complete_strains.csv
write.csv(annotated, file = "~/GitHub/riboSeed/sample_data/all_complete_strains.csv", row.names = F, quote = F)
'
write("USAGE: Rscript get_closest_strain.R path/to/full_strains.csv genus species ramdom_seed(optional)", stderr())
# data from ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
args <- commandArgs(trailingOnly = TRUE)
combined_complete_genomes <- as.character(args[1])
# combined_complete_genomes <- "~/GitHub/riboSeed/sample_data/all_complete_strains.csv"

annotated <- read.csv2(combined_complete_genomes, header = T, stringsAsFactors = F, sep=",")
# 

genus <-  as.character(args[2])
species <- as.character(args[3])
if (species == "") stop("must provide species")
random <- as.numeric(args[4])
# random = "-r"
# genus = "Rhodobacter"
# species = "Sphaeroides"
RANDOM = !is.na(random)

write(random, stderr())

make_subset <- function(genus, species, annotated){
  # check if genus hits
  genus_subset <- annotated[grepl(genus, annotated$genus, ignore.case = T),]
  if (nrow(genus_subset) == 1){  # return if single genus hit
    return(genus_subset)
  } else if (nrow(genus_subset) == 0){ # if no genus hits, stop
    stop("No hits at the Genus level")
  } else { # if muti=ple hits, subset and start looking at species level
    species_subset <- genus_subset[grepl(species, genus_subset$species, ignore.case = T),]
    if (nrow(species_subset) == 1){  # if sinlge species hit, return it
      return(species_subset)
    } else if (nrow(species_subset) == 0){ # if no species hits, return random genus levels
      stop("No hits at the species level")
    } else { # ie, there are multiple species hits, check for repre
      species_repr_subset <- species_subset[species_subset$Reference == "REFR",]
      if (nrow(species_repr_subset) == 1){
        return(species_repr_subset)
      } else if (nrow(species_repr_subset) > 1){
        return(species_repr_subset)
      } else { # if no representative strains, either return alla random one
        return(species_subset)
      }
    }
  }
}
if (RANDOM) set.seed(random)
subset_genus_species <- make_subset(genus, species, annotated)
string <- subset_genus_species$Replicons
get_chromosomes = function(string){
  chroms = unlist(strsplit(x = string, split="chromosome"))
  accs = gsub(":(.*?)/(.*)", "\\1", chroms)
  return(accs[accs != ""])
}
write(get_chromosomes(string), stdout())

write.csv("subset_of_complete_strains.csv", x = subset_genus_species)
