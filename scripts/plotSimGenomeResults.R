# version 0.0.2
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("msa")
# biocLite("ggtree")
library(ggplot2)
library(reshape2)
library(ape)
library(ggtree)
require(msa)
library(seqinr)
require(phytools)
require(ggstance)
analysis_date = "2017-06-15"
src_folders <-"~/GitHub/riboSeed/manuscript_results/simulated_genome/"

names = c("sourcepath", "assembly", "total", "good",  "ambiguous", "bad")
dir(src_folders)
folders<- dir(src_folders, pattern = "simulatedGenomeResults_",full.names = T)
folders
folders = folders[1:8]

out_folder <- paste0("~/GitHub/riboSeed/", analysis_date, "-simulated-genome/")
dir.create(out_folder)
results_table <- as.data.frame(setNames(replicate(6,numeric(0), simplify = F), names))
for (folder in folders){
  try(
    results_table <- rbind(results_table,
                           read.csv2(
                             paste0(folder, "/riboScore/riboScore_report.txt"), 
                             sep="\t", stringsAsFactors = F, header = F, col.names = names))
  )
}
str(results_table)
results_table$name <-gsub("\\.fasta", "", results_table$assembly)
results_table$assembly <- NULL
  
print(results_table)


tall <- melt(results_table, id.vars = c("sourcepath", "name", "total"), measure.vars = c("good", "bad", "ambiguous"), factorsAsStrings = T)
str(tall)
tall$name <- as.factor(tall$name)
levels(tall$name) <-c("E. Coli", "K. pneumoneae")
tall$prettynames <- levels(tall$name) 
boxy <- ggplot(tall, aes(x=variable, y=value)) + facet_grid(~prettynames) + 
  # geom_vline(0, xintercept = 0, color="grey30") + 
  geom_boxplot(outlier.colour =  NA,color="grey40", fill="grey60",width=.7)+
  geom_point(position=position_jitter(width = .02, height = .12), shape=1)+
  scale_x_discrete(expand=c(.05,.05))+
  scale_y_continuous(expand=c(.01,.01), limits = c(-.2,7.5), breaks = 0:7)+
  theme(
        # axis.title.y = element_blank(),
        # axis.text.y = element_blank(), 
        # axis.ticks = element_blank(), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        
        plot.title = element_text(size =15))+
  labs(y="Count", x="rDNA Assembly ",  title="")

boxy
plot1 <- paste0(out_folder, "plot1")

pdf(file = paste0(plot1, ".pdf"), width = 4, height = 4)
boxy
dev.off()

png(filename = paste0(plot1, ".png"), width = 4, height = 4, units = "in", res = 300)
boxy
dev.off()


