# version 0.0.3
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("msa")
# biocLite("ggtree")
library(ggplot2)
library(reshape2)
library(dplyr)

help = "
This takes the path to a directory containing two subdirs of files: the hits
files of coli, and of kleb.  These files should not be cat-ed into one, as it
is easier to fill iun missing values if they are kept separate.
"

src_folder <-"/home/nicholas/GitHub/riboSeed/2017-10-06-coli_kleb_simulated_reads"
names = c("sourcepath", "variable", "start",  "end1", "end2", "localpath")
dir(src_folder)
files <- dir(src_folder, pattern = "*.txt",full.names = T, recursive = T)
files

out_folder <- paste0("~/GitHub/riboSeed/2017-10-09-simReads/")

dir.create(out_folder)
results_table <- as.data.frame(
  setNames(replicate(7, numeric(0), simplify = F), c(names, "num")))
for (file in files){
  # read in the file
  new <- read.csv2(
    file, sep="\t", stringsAsFactors = F, header = F, col.names = names)
  # this is kind of a sillly variable that gets used for counting! is 1 if   
  # data exists, and 0 later if we fill in for missing data
  new$num <- 1
  # fill in missing values with 0's
  for (var in c("good","bad", "?")){
    if ( var %in% new$variable){
      print("found")
    } else{
      print(paste("missing", var))
      new <- rbind(
        new, 
        data.frame("localpath"=new$localpath,
                   "sourcepath" = new$sourcepath[1], 
                   "variable" = var, 
                   "start" = "NotApplicable",
                   "end1" = "NotApplicable",
                   "end2" = "NotApplicable",
                   "num" = 0)
      )        
    }
  }
  new$localpath <- file
  results_table <- rbind(results_table, new)
}

str(results_table)
# parse replicate number
results_table$rep <- as.numeric(gsub(".*rep(.).*", "\\1", results_table$sourcepath))
# parse coverage depth
results_table$coverage <- as.numeric(gsub(".*COV_(.*?)_.*", "\\1", results_table$sourcepath))

# parse bug
results_table$bug <- ifelse(grepl("/coli/", results_table$localpath), "E. coli", "K. pneumoneae")

# parse technology
results_table$tech <- (gsub(".*sim_(.*?)_COV.*", "\\1", results_table$sourcepath))

results_table$sourcepath <-NULL
results_table$localpath <-NULL
#results_table$name <-gsub("\\.fasta", "", results_table$sourcepath)
tall <- results_table[, !colnames(results_table) %in% c("start", "end1", "end2")] %>% 
  group_by(coverage, tech, variable, rep, bug) %>%
  mutate(count=sum(num)) %>%
  as.data.frame()
tall <- tall[!duplicated(tall),colnames(tall) != "num" ]
View(tall)

ggplot(tall, aes(x=variable, y=count)) + 
  geom_boxplot() + geom_point(position=position_jitter(width=.5, height=0)) +
  facet_grid(~tech+coverage+bug)

tall$val = tall$count
# set things we need to rename as factors
tall$variable <- as.factor(tall$variable)
tall$tech<-as.factor(tall$tech)
tall$coverage<-as.factor(tall$coverage)
# set the factor levels to the proper labels
levels(tall$tech)<-c("GenomeAnalyzerII, 50bp", "GenomeAnalyzerII, 75bp", "HiSeq, 100bp")
levels(tall$coverage)<-c("20x Coverage", "50x Coverage")
levels(tall$variable) <- c("Ambiguous", "Incorrect", "Correct")


# plot theme
theme_sim <- function(){
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
      axis.line = element_line(color = 'black', size = .1)
    )
    
}
(coliplot<- ggplot(tall[tall$bug=="E. coli",], 
                   aes(x=variable, y=val)) + 
    geom_boxplot(outlier.colour =  NA,color="grey20", fill="grey65",width=.7)+
    geom_point(position=position_jitter(width=.1, height=.1), shape=1) +
    facet_grid(coverage~tech)+
    scale_x_discrete(expand=c(.05,.05), breaks=c("Correct", "Incorrect", "Ambiguous"))+
    scale_y_continuous(expand=c(.01,.01), limits = c(-.2,7.5), breaks = 0:7)+
    theme_sim() +
    labs(y="Count", x=" ",  title=""))
set.seed(27)
(klebplot<- ggplot(tall[tall$bug=="K. pneumoneae",], aes(x=variable, y=val)) + 
    geom_boxplot(outlier.colour =  NA,color="grey20", fill="grey65",width=.7)+
    geom_point(position=position_jitter(width=.1, height=.10), shape=1) +
    facet_grid(coverage~tech)+
    scale_x_discrete(expand=c(.05,.05))+
    scale_y_continuous(expand=c(.01,.01), limits = c(-.2,8.5), breaks = 0:8)+
    theme_sim()+
    labs(y="Count", x=" ",  title="")
  )
  
(combplot2<- ggplot(tall, aes(x=tech, y=val, fill=variable, color=variable)) + 
    geom_boxplot(outlier.colour =  NA, width=.7, alpha=0.6)+
    geom_point(position=position_jitterdodge(jitter.width=.1, jitter.height=.1), shape=1, color="black") +
    facet_grid(bug~coverage)+#, scales = "free_y")+
    scale_x_discrete(expand=c(.05,.05))+
    scale_color_discrete(guide="none")+
    scale_y_continuous(expand=c(.01,.01), limits = c(0,8), breaks = 0:8)+
    theme_bw() +
    #eliminates background, gridlines, and chart border
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color="black", size = .1),
      strip.background = element_rect(colour = "grey85", size = 2),
      strip.text = element_text(size = 9, face="italic"),
      axis.text.x = element_text(angle=45, hjust=1, size = 12),
      axis.text.y = element_text(angle=0, size = 12),
      axis.title  = element_text(angle=0, size = 13),
      axis.line = element_line(color = 'black', size = .1)) +
    labs(y="Count", x="Technology, Read Length ",  title="", fill="rDNA Assembly"))

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




