#!/usr/bin/env Rscript
# version 0.0.2
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("msa")
# biocLite("ggtree")
#  Imports are moved til after argument parses!
require(argparse, quietly = T, warn.conflicts = F)
help <- "
This script is used to plot the results from runDegenerate.  As an input, it takes a file containing the combined ribo score reports, and requires an output directory.

USAGE:
Rscript plotDegenResults.R -r path/to/combinedresults.txt -o ./output_dir/

"
###################################################################
# create parser object
parser <- ArgumentParser()
# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-r", "--combined_reports", action="store",
                    dest="combined_reports", help="")
# parser$add_argument("-r", "--results_csv", action="store", 
#                     default="NA",
#                     dest="results", help="columns must be name, seed, value, variable")
parser$add_argument("-o", "--out_dir", action="store",
                    dest="out_folder", help="output_directory")

args <- parser$parse_args()
if (is.null(args$out_folder)){
  message(help)
  stop("Must provide output path")
}
print(args)
if (dir.exists(args$out_folder)){
  print("output dir exists!")
  quit()
}
require(ggplot2)
require(reshape2, quietly = T, warn.conflicts = F)
require(dplyr)
require(tidyr)

# (analysis_date <- gsub("(.*)([0-9]{4}-[0-9]{2}-[0-9]{2})(.*)", "\\2", args$src_folder))

(names <- c("sourcepath", "assembly", "total", "good",  "ambiguous", "bad"))

dir.create(args$out_folder)
reportsdf <- as.data.frame(setNames(replicate(6,numeric(0), simplify = F), names))

freqs = c(0.0, 0.0025, 0.0050, 0.0075, 0.0100, 0.0150, 0.0200, 0.0250, 0.0500, 0.0750, 0.1000, 0.1250, 0.1500, 0.1750, 0.2000, 0.2250, 0.2500, 0.2750, 0.3000)
str(reportsdf)
temp_new_table <- read.csv2(
  args$combined_reports, sep="\t", stringsAsFactors = F, header = F, col.names = names)
# temp_new_table <- read.csv2("~/GitHub/riboSeed/2017-10-degen-results/2017-10-13-combined_reports.txt", sep="\t", stringsAsFactors = F, header = F, col.names = names)

(temp_new_table$freq <- as.numeric(gsub("(.*)/seed_(.*?)/mauve", "\\2", temp_new_table$sourcepath)))
(temp_new_table$where <- factor(gsub("(.*)_degenerate_output_(.*?)_(.*)", "\\2", temp_new_table$sourcepath)))
(temp_new_table$seed <- as.numeric(gsub("(.*)_degenerate_output_(.*?)_(.*?)/(.*)", "\\3", temp_new_table$sourcepath)))

# remove de novo assemblies!
temp_new_table <- temp_new_table[grepl("de_fere_novo", temp_new_table$assembly),]

temp_new_table$sourcepath <-NULL
temp_new_table$assembly <-NULL

seeds <- unique(temp_new_table$seed)
print(seeds)
print(freqs)
str(temp_new_table)
table(temp_new_table$where)

for (where in c("ALL", "FLANK")){
  for (seed in seeds){
    subdf <- temp_new_table[temp_new_table$seed == seed &
                            temp_new_table$where == where, ]
    missings <- freqs[!(freqs %in% subdf$freq)]
    if (!length(missings) == 0){
      print(paste(where, "and", seed, "missing", paste(missings, collapse=" ")))
      dummydf <-data.frame(
        "total" = max(subdf$total),
        "good" = 0,
        "ambiguous" = 0,
        "bad" = 0,
        "freq" = missings,
        "seed" = seed,
        "where" = where)
      temp_new_table <- rbind(temp_new_table, dummydf)
    }
  }
}
  
reportsdf <- temp_new_table
write.csv(reportsdf, file = file.path(args$out_folder, "test"))
str(reportsdf)
reportsdf$value <- reportsdf$good
tall <- reportsdf[,!colnames(reportsdf) %in% c("bad",  "ambiguous", "total", "good") ]

levels(tall$where) <- c("rDNA + Flanking", "Flanking")
#  from https://stackoverflow.com/questions/3010403/jitter-if-multiple-outliers-in-ggplot2-boxplot
require(dplyr)
id_outliers <- function(x){
  q <- quantile(x,c(0.25,0.75))
  iqr <- abs(diff(q))
  ifelse((x < q[1] - 1.5*iqr) | (x > q[2] + 1.5*iqr), T, F)
}
tall <- tall %>%
  group_by(where, freq) %>%
  mutate(out= id_outliers(value)) %>% 
  as.data.frame()

sumary_tall <- tall %>%
  group_by(where, freq) %>%
  count(value)

str(tall)

labelsdf <- data.frame(x1=c(0,.05, .15), 
                       x2=c(.05, .15, .31), 
                       y1=c(0,0,0),
                       y2=c(7,7,7), labs=c("Species", "Genus","Other" ))
(line_log <- ggplot(tall[tall$freq != 0,],
                    aes(x=freq, color=where, y=value))+ 
    annotate("rect", xmin=min(tall[tall$freq != 0, "freq"])*.9, xmax=0.05, ymin=0, ymax=8, alpha=0.2, fill="green", color=NA)+
    annotate("text", x=0.008,  y=7.5, label="Intrapecies")+
    annotate("text", x=0.09,  y=7.5, label="Interspecies")+
    geom_smooth(aes(group=where, color=where)) +
    geom_point(data=tall[tall$out==T, ], color="black", alpha=.2, position=position_jitter(width = .01, height = .05))+
    geom_point(data=tall[tall$out==F, ],  alpha=.1, position=position_jitter(width = .03, height = .05), aes(x=freq, y=value))+
    geom_hline(0, yintercept = 0, color="grey20") + 
    scale_x_log10(limits = c(NA, 0.3), expand=c(.01, .01), breaks=c(.01, .02, .05, .1, .15, .2, .3))+
    scale_y_continuous(expand=c(.03, .03), limits = c(-.05, 8), breaks = 0:7)+
    scale_color_manual(values=c("red", "blue"))+
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill=NA),
      strip.text = element_text( angle = 0, face="bold", size=10 ),
      plot.title = element_text(size=15))+
    labs(y="Correctly-assembled rDNAs", x="Substitution Frequency", title="", color="Permitted Substitutions")
)
labelss <- c()
for (freq in freqs){
  if (freq %in% c(0.0, .05, 0.1, 0.2, 0.3)){
    labelss <- c(labelss, freq)
  } else{
    labelss <- c(labelss, "")
  }
}

line_lin <- ggplot(tall,#[tall$freq != 0, ],
                    aes(x=freq, color=where, y=value))+ 
    annotate("rect", xmin=min(tall$freq), xmax=0.05, ymin=0, ymax=8, 
             alpha=0.15, fill="purple", color=NA)+
    annotate("text", x=0.025,  y=7.5, label="Intraspecies", size=5)+
    annotate("text", x=0.090,  y=7.5, label="Interspecies", size=5)+
    geom_smooth(size=0.5,
        method="loess", 
      formula = y ~ x, span=1)+
    geom_point(data=tall[tall$out==T, ], shape=3, size=1, alpha=.8, position=position_jitter(width = .0005, height = .05))+
    geom_point(data=tall[tall$out==F, ], shape=1, size=1, alpha=.5, position=position_jitter(width = .001, height = .05))+
    geom_hline(0, yintercept = 0, color="grey20") + 
    scale_x_continuous(limits = c(0, 0.3), expand=c(0.01, 0), breaks=freqs, labels =labelss)+
    scale_y_continuous(expand=c(0, 0), limits = c(-0.4, 8), breaks = 0:7)+
    scale_color_manual(values=c("#E69016", "#2A61AD"))+
    theme(
      axis.title = element_text( angle = 0, size=17 ),
      plot.background = element_blank(),
      panel.background = element_blank(),
      legend.position = "bottom",
      # legend.key = element_rect(fill = "transparent", color="transpar"),
      strip.background = element_rect(fill=NA),
      strip.text = element_text( angle = 0, face="bold", size=10 ),
      plot.title = element_text(size=15))+
    labs(y="Correctly-assembled rDNAs", x="Substitution Frequency", title="", 
         # color="Permitted Substitutions")
          color="")
## oh joy...
library(ggridges)
ggplot(tall, aes(x=value, y=factor(freq, levels=unique(freq)), color=where, fill=where)) + geom_density_ridges(alpha=.5) + scale_y_discrete(limits = as.character(rev(unique(tall$freq))))

###################################################
# For summarized table
###################################################

sum_line_lin <- ggplot(tall,#[tall$freq != 0, ],
                   aes(x=freq, color=where, y=value))+ 
  annotate("rect", xmin=min(tall$freq), xmax=0.05, ymin=0, ymax=8, 
           alpha=0.15, fill="purple", color=NA)+
  annotate("text", x=0.025,  y=7.5, label="Intraspecies", size=5)+
  annotate("text", x=0.090,  y=7.5, label="Interspecies", size=5)+
  geom_smooth(size=0.5,
              method="loess", 
              formula = y ~ x, span=1)+
  geom_point(data=sumary_tall, aes(size=n, fill=where), shape=21, alpha=.4)+#, position=position_jitter(width = .001, height = .05))+
  geom_hline(0, yintercept = 0, color="grey20") + 
  scale_x_continuous(limits = c(0, 0.3), expand=c(0.02, 0), breaks=freqs, labels =labelss)+
  scale_y_continuous(expand=c(0, 0), limits = c(-0.4, 8), breaks = 0:7)+
  scale_color_manual(values=c("#FF8000", "#0048AD"))+
  scale_fill_manual(values=c("#FF8000", "#0048AD"))+
#  scale_size_continuous(range=c(.25,8)) +
  scale_size_area(max_size=8, breaks=c(1,25,50,100)) +
  theme(
    axis.title = element_text( angle = 0, size=17 ),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom",
    legend.key = element_blank(),
    # legend.key = element_rect(fill = "transparent", color="transpar"),
    strip.background = element_rect(fill=NA),
    strip.text = element_text( angle = 0, face="bold", size=10 ),
    plot.title = element_text(size=15))+
  labs(y="Correctly-assembled rDNAs", x="Substitution Frequency", title="",
       size="", 
       # color="Permitted Substitutions")
       color="", fill="")

###################################################


#  coord_flip()
writePlot <- function(pl, file, h, w){
  pdf(file = paste0(file, ".pdf"), width = w, height = h)
  print(pl)
  dev.off()
  png(filename = paste0(file, ".png"), width = w, height = h, units = "in", res = 300)
  print(pl)
  dev.off()
}


writePlot(
  pl=line_lin, 
  file=paste0(args$out_folder, "degenerate_lineplot"),
  w=8, h=6)
writePlot(
  pl=sum_line_lin, 
  file=paste0(args$out_folder, "degenerate_bubble_lineplot"),
  w=8, h=6)


