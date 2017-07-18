# version 0.0.2
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("msa")
# biocLite("ggtree")
require(argparse, quietly = T, warn.conflicts = F)
library(ggplot2)
library(reshape2)
library(dplyr)




combined_reports_kleb <-"/home/nicholas/GitHub/riboSeed/manuscript_results/simulated_reads_kleb/kleb_combined_reports.txt"
combined_reports_coli <-"/home/nicholas/GitHub/riboSeed/manuscript_results/simulated_reads_coli/2017-06-23/coli_proper_combined_reports.txt"
(names <- c("sourcepath", "assembly", "total", "Correct",  "Ambiguous", "Incorrect"))

out_folder <- paste0("~/GitHub/riboSeed/2017-06-23-simReads_replot/")

dir.create(out_folder)
results_table <- as.data.frame(
  setNames(replicate(7, numeric(0), simplify = F), c(names, "bug")))
colidf <- read.csv2(
  combined_reports_coli, sep="\t", stringsAsFactors = F, header = F, col.names = names)
colidf$bug<-"E. coli"
klebdf <- read.csv2(
  combined_reports_kleb, sep="\t", stringsAsFactors = F, header = F, col.names = names)
klebdf$bug<-"K. pneumoneae"
results_table <- rbind(results_table, colidf, klebdf)

str(results_table)
results_table$rep <- as.numeric(gsub(".*rep(.).*", "\\1", results_table$sourcepath))
results_table$coverage <- as.numeric(gsub(".*COV_(.*?)_.*", "\\1", results_table$sourcepath))

results_table$tech <- (gsub(".*sim_(.*?)_COV.*", "\\1", results_table$sourcepath))
results_table$sourcepath <-NULL
str(results_table)

tall <- melt(
  data = results_table[results_table$assembly == "riboSeed_de_fere_novo_contigs.fasta",],
  id.vars = c("coverage", "tech", "rep", "bug"), measure.vars = c("Correct", "Incorrect", "Ambiguous"),
  variable.name = "variable", value.name = "val")
tall$tech<-as.factor(tall$tech)
tall$coverage<-as.factor(tall$coverage)
levels(tall$tech)<-c("GenomeAnalyzerII, 50bp", "GenomeAnalyzerII, 75bp", "HiSeq, 100bp")
levels(tall$coverage)<-c("20x Coverage", "50x Coverage")
(coliplot<- ggplot(tall[tall$bug=="E. coli",], aes(x=variable, y=val)) + 
    geom_boxplot(outlier.colour =  NA,color="grey40", fill="grey60",width=.7)+
    geom_point(position=position_jitter(width=.1, height=.1), shape=1) +
    facet_grid(coverage~tech)+
    scale_x_discrete(expand=c(.05,.05))+
    scale_y_continuous(expand=c(.01,.01), limits = c(-.2,7.5), breaks = 0:7)+
    theme_bw() +
    #eliminates background, gridlines, and chart border
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color="black", size = .1),
      strip.background = element_rect(colour = "grey85", size = 2),
      strip.text = element_text(size = 9),
      axis.text.x = element_text(angle=45, hjust=1, size = 12),
      axis.text.y = element_text(angle=0, size = 12),
      axis.title  = element_text(angle=0, size = 13),
      axis.line = element_line(color = 'black', size = .1)) +
    labs(y="Count", x=" ",  title=""))
set.seed(27)
(combplot2<- ggplot(tall, aes(x=tech, y=val, fill=variable, color=variable)) + 
    geom_boxplot(outlier.colour =  NA, width=.7, alpha=0.5)+
    geom_point(position=position_jitterdodge(jitter.width=.1, jitter.height=.1), shape=1, color="black") +
    facet_grid(bug~coverage)+
    scale_x_discrete(expand=c(.05,.05))+
    scale_color_discrete(guide="none")+
    scale_y_continuous(expand=c(.01,.01), limits = c(-.2,7.5), breaks = 0:7)+
    theme_bw() +
    #eliminates background, gridlines, and chart border
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color="black", size = .1),
      strip.background = element_rect(colour = "grey85", size = 2),
      strip.text = element_text(size = 9),
      axis.text.x = element_text(angle=45, hjust=1, size = 12),
      axis.text.y = element_text(angle=0, size = 12),
      axis.title  = element_text(angle=0, size = 13),
      axis.line = element_line(color = 'black', size = .1)) +
    labs(y="Count", x="Technology, Read Length ",  title="", fill="rDNA Assembly"))

(klebplot<- ggplot(tall[tall$bug=="K. pneumoneae",], aes(x=variable, y=val)) + 
  geom_boxplot(outlier.colour =  NA,color="grey40", fill="grey60",width=.7)+
  geom_point(position=position_jitter(width=.1, height=.10), shape=1) +
  facet_grid(coverage~tech)+
  scale_x_discrete(expand=c(.05,.05))+
  scale_y_continuous(expand=c(.01,.01), limits = c(-.2,8.5), breaks = 0:8)+
    theme_bw() +
    #eliminates background, gridlines, and chart border
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color="black", size = .1),
      strip.background = element_rect(colour = "grey85", size = 2),
      strip.text = element_text(size = 9),
      axis.text.x = element_text(angle=45, hjust=1, size = 12),
      axis.text.y = element_text(angle=0, size = 12),
      axis.title  = element_text(angle=0, size = 13),
      axis.line = element_line(color = 'black', size = .1)) +
    labs(y="Count", x=" ",  title=""))

stattab <- tall %>%
  group_by(coverage, variable, tech, bug) %>%
  mutate(mean = mean(val), stdev=sd(val)/sqrt(5)) %>%
  as.data.frame()

options(scipen = 2)
(stattab$textout <- paste0(stattab$mean, " (", stattab$stdev, ")"))
stattab <- stattab[, !colnames(stattab) %in% c("val", "rep", "mean", "stdev")]
stattab$coverage <-as.factor(stattab$coverage)
stattab$tech <-as.factor(stattab$tech)
stattab$variable <-as.factor(stattab$variable)
stattab <- stattab[!duplicated(stattab),]
stattab
# dcast(stattab, coverage + variable ~ tech, value.var = "val")
# str(stattab)

require(tidyr)
stab <- spread(stattab, key = variable, value = textout)
stab

write.csv(stab, file = file.path(out_folder,"combined.csv"))


pdf(file = paste0(out_folder, "coli.pdf"), width = 5.5, height = 6.5)
coliplot
dev.off()
pdf(file = paste0(out_folder, "kleb.pdf"), width = 5.5, height = 6.5)
klebplot
dev.off()
pdf(file = paste0(out_folder, "simulated_reads_combined.pdf"), width = 6, height = 6.5)
combplot2
dev.off()
png(filename = paste0(out_folder, "kleb.png"), width = 5.5, height = 6.5, units = "in", res = 300)
klebplot
dev.off()
png(filename = paste0(out_folder, "coli.png"), width = 5.5, height = 6.5, units = "in", res = 300)
coliplot
dev.off()


