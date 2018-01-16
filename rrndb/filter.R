library(dplyr)
file = "~/rrnDB-5.4.tsv"
rrndb <- read.csv(file, sep="\t", stringsAsFactors = F)
rrndb$Data.source.lineage <- gsub("(.*)RDP_ID; (.*)", "\\2", rrndb$Data.source.lineage)
rrndb$Data.source.lineage <-ifelse(
  "" == rrndb$Data.source.lineage, 
  rrndb$RDP.taxonomic.lineage,
  rrndb$Data.source.lineage
)
i=1
goodnames = c("domain", "phylum", "order", "family", "genus")
new_df <- data.frame(domain=NA, phylum=NA, order=NA, family=NA, genus=NA) 
for (i in 1:nrow(rrndb)){
  fields = strsplit(rrndb[i, "Data.source.lineage"], ";", fixed = T)[[1]]
  names(fields) <- gsub(" ", "", gsub("(.*)\\|(.*)", "\\2", fields))
  fields <- gsub(" ", "", gsub("(.*)\\|(.*)", "\\1", fields))
  fields <- gsub('\"', "", fields, fixed = T)
  new_df = rbind(new_df, fields[goodnames])
}
new_df <- new_df[2:nrow(new_df), ]
new_df$strain <- rrndb$Data.source.organism.name
# sue me
new_df$species <- gsub("(.*?)\\s(.*?)\\s.*", "\\2", paste0(new_df$strain, " "))
new_df$bioproject <- rrndb$BioProject
new_df$rdnas <- ifelse(is.na(rrndb$X23S.gene.count),
                       rrndb$X16S.gene.count,
                       rrndb$X23S.gene.count)

summary <- new_df %>%
  group_by(phylum) %>%
  filter(max(rdnas) > 1) %>%
  filter(n()>1) %>%
  summarise(count = n(), max=max(rdnas), min=min(rdnas))

summary_full <- new_df %>%
  filter(rdnas > 1) %>%
  filter(!is.na(genus)) %>%
  filter(!is.na(bioproject)) %>%
  group_by(phylum, order, family, genus, species) %>%
  filter(n()>2) %>%
  summarise(count = n(), max=max(rdnas), min=min(rdnas)) 

summary <- new_df %>%
  filter(rdnas > 1) %>%
  filter(!is.na(genus)) %>%
  filter(!is.na(bioproject)) %>%
  group_by(phylum, order, family, genus, species) %>%
  filter(n()>2) %>%
  as.data.frame() %>%
  group_by(phylum) %>%
  summarise(count = n(), max=max(rdnas), min=min(rdnas)) 

