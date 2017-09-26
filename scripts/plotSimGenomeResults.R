# version 0.0.2
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("msa")
# biocLite("ggtree")
library(ggplot2)
library(reshape2)

analysis_date <- Sys.Date()
src_folders <-"~/GitHub/riboSeed/manuscript_results/simulated_genome/"

names = c("sourcepath", "assembly", "total", "Correct",  "Ambiguous", "Incorrect")
dir(src_folders)
folders<- dir(src_folders, pattern = "simulatedGenomeResults_",full.names = T)
folders
# folders = folders[1:5]

out_folder <- paste0("~/GitHub/riboSeed/", analysis_date, "-simulated-genome-replot/")
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
#  add in kleb de fere novo values;  
kleb_de_fere_novo <- data.frame(
  "sourcepath"=paste0("manually_added", 1:10),
  "assembly"=rep("kleb_de_fere_novo.fasta", 10),
  "total"=rep(7, 10), 
  "Correct"=rep(0, 10),
  "Ambiguous"=rep(0, 10),
  "Incorrect"=rep(0, 10))

results_table <- rbind(results_table, kleb_de_fere_novo)
str(results_table)
results_table$name <-gsub("\\.fasta", "", results_table$assembly)
results_table$assembly <- NULL
print(results_table)

tall <- melt(results_table, id.vars = c("sourcepath", "name", "total"), measure.vars = c("Correct", "Incorrect", "Ambiguous"), factorsAsStrings = F)
str(tall)
tall$name <- factor(tall$name,  levels=c("coli_de_fere_novo", "coli_de_novo", "kleb_de_fere_novo"))
levels(tall$name) <-c("E. coli (de fere novo)", "E. coli (de novo)", "K. pneumoniae (de fere novo)")
tall$name <- factor(tall$name, levels=levels(tall$name)[c(2,1,3)]) 

(boxy <- ggplot(tall, aes(x=variable, y=value)) + facet_grid(~name) + 
  # geom_vline(0, xintercept = 0, color="grey30") + 
  geom_boxplot(outlier.colour =  NA,color="grey40", fill="grey60",width=.7)+
  geom_point(position=position_jitter(width = .02, height = .12), shape=1)+
  scale_x_discrete(expand=c(.05,.05))+
  scale_y_continuous(expand=c(.01,.01), limits = c(-.2,7.5), breaks = 0:7)+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", size = .1),
    strip.background = element_rect(colour = NA),
    strip.text.x =element_text(face="italic"),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.line = element_line(color = 'black', size = .1)) +
  labs(y="Count", x="rDNA Assembly ",  title=""))

boxy

boxy_group_strain <- ggplot(tall, aes(x=name, y=value, color=variable, fill=variable)) +
  # geom_vline(0, xintercept = 0, color="grey30") + 
  geom_boxplot(outlier.colour =  NA, width=.7, alpha=0.5)+
  geom_point(position=position_jitterdodge(dodge.width = .7, jitter.width=.1, jitter.height=.1 ), shape=1, color="black")+
  scale_color_discrete(guide="none")+
  scale_x_discrete(expand=c(.05,.05))+
  scale_y_continuous(expand=c(.01,.01), limits = c(-.2,7.5), breaks = 0:7)+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", size = .1),
    strip.background = element_rect(colour = NA),
    axis.text.x = element_text(angle=45, hjust=1),
#    legend.position=c(0.7, .7),
    axis.line = element_line(color = 'black', size = .1)) +
  labs(y="Count", x="Strain (method)",  title="", fill="rDNA assembly")

boxy_group_strain
boxy_group_rdna <- ggplot(tall, aes(x=variable, y=value, fill=name, color=name)) +
  # geom_vline(0, xintercept = 0, color="grey30") + 
  geom_boxplot(outlier.colour =  NA, width=.7, alpha=0.5)+
  geom_point(position=position_jitterdodge(dodge.width = .7, jitter.width=.1, jitter.height=.1 ), shape=1, color="black")+
  scale_fill_discrete(position="bottom")+
  scale_color_discrete(guide="none")+
  scale_x_discrete(expand=c(.05,.05))+
  scale_y_continuous(expand=c(.01,.01), limits = c(-.2,7.5), breaks = 0:7)+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", size = .1),
    strip.background = element_rect(colour = NA),
    axis.text.x = element_text(angle=45, hjust=1),
#    legend.position=c(0.7, .7),
    axis.line = element_line(color = 'black', size = .1)) +
  labs(y="Count", x="rDNA Assembly ",  title="", fill="Strain (method)")

boxy_group_rdna

pdf(file = file.path(out_folder, "simulated_genome.pdf"), width = 6, height = 5)
boxy
# boxy_group_rdna
# boxy_group_strain
dev.off()
png(filename = paste0(out_folder, "simulated_genome.png"), width = 4, height = 4, units = "in", res = 300)
boxy
dev.off()


