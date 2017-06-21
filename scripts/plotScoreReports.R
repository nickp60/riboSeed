# version 0.0.2
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("msa")
# biocLite("ggtree")
require(argparse, quietly = T, warn.conflicts = F)
library(ggplot2)
library(reshape2)
library(dplyr)




combined_reports <-"/home/nicholas/GitHub/riboSeed/manuscript_results/simulated_reads_kleb/kleb_combined_reports.txt"
combined_reports <-"/home/nicholas/GitHub/riboSeed/manuscript_results/simulated_reads_coli/coli_combined_reports.txt"
(names <- c("sourcepath", "assembly", "total", "good",  "ambiguous", "bad"))

out_folder <- paste0("~/GitHub/riboSeed/2017-06-16-simReads/")

dir.create(out_folder)
results_table <- as.data.frame(
  setNames(replicate(6, numeric(0), simplify = F), c(names)))
new <- read.csv2(
    combined_reports, sep="\t", stringsAsFactors = F, header = F, col.names = names)
results_table <- rbind(results_table, new)

str(results_table)
results_table$rep <- as.numeric(gsub(".*rep(.).*", "\\1", results_table$sourcepath))
results_table$coverage <- as.numeric(gsub(".*COV_(.*?)_.*", "\\1", results_table$sourcepath))

results_table$tech <- (gsub(".*sim_(.*?)_COV.*", "\\1", results_table$sourcepath))
results_table$sourcepath <-NULL
str(results_table)
#results_table$name <-gsub("\\.fasta", "", results_table$sourcepath)
#tall <- results_table[, !colnames(results_table) %in% c("start", "end1", "end2")] %>%
tall <- melt(data = results_table[results_table$assembly == "riboSeed_de_fere_novo_contigs.fasta",], id.vars = c("coverage", "tech", "rep"), measure.vars = c("good", "bad", "ambiguous"), variable.name = "variable", value.name = "val")

ggplot(tall, aes(x=variable, y=val, fill=as.factor(coverage))) + 
  geom_boxplot() + #geom_point()+#position=position_jitter(width=.1, height=0)) +
  facet_grid(~tech)

stattab <- tall %>%
  group_by(coverage, variable, tech) %>%
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

ftable(stattab$variable, stattab$coverage, stattab$tech, 
       row.vars = c(1,2), 
       col.vars = c(3), stattab$textout)
require(tidyr)
require(xtable)
stab <- spread(stattab, key = variable, value = textout)
stab
(outf  <- gsub(paste0(".*",file.path("", ""), "(.*?)_combined.*"), "\\1", combined_reports))
write.csv(stab, file = file.path(out_folder, paste0(outf, ".csv")))
xftbl <- xtable(stab) 
print.xtableFtable(xftbl)


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


