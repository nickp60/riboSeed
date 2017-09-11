#!/usr/bin/env Rscript
#############################################################################################
help<-"
Blast  Parser
20170103
"
print("USAGE: path/to/snag/output/ path/to/output/dir/ threshold_percentage")
require(dplyr)
#require(reshape2)
require(ggplot2)
this_version<-'0.3'
if (interactive()){
  args= c("~/GitHub/riboSeed/manuscript_results/entropy/sakai_snag_mafft",
          "~/GitHub/riboSeed/BLAST_results/", 90)
}else{
  args <- commandArgs(trailingOnly = TRUE)
}
snagdir <- as.character(args[1])
outdir <- as.character(args[2])
threshold <- as.character(args[3])
VERBOSE_bool<-T

blast_results <- file.path(snagdir, "BLAST_results",
  grep("*_results_merged.csv", 
       dir(file.path(snagdir,"BLAST_results")), value=T))

datetime<-gsub("[:-]","",gsub("([\\d:-]*)(\\D*)","\\1",Sys.time()))

dir.create(outdir)
##########

#read in csv
data <- read.csv2(
  blast_results, header=F, stringsAsFactors = F, sep="\t", 
  col.names = c(
    "query_id", "subject_id", "identity_perc", "alignment_length", 
    "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end",
    "evalue", "bit_score"))
data$flank <-as.numeric(gsub("(.*)_(\\d*)-bp_flanks", "\\2", data$query_id))
data$region <-(gsub("(.*?)_(.*?)_(.*?)_.*", "\\3", data$query_id))

print("DATA:")
print(str(data))
print("FLANKING REGIONS")
print(unique(data$flank))
print("REGIONS")
print(unique(data$region))


data2<-data %>%
  # for each query_id, or region
  group_by(query_id) %>%
  # normalize by length of longest match, which should be the "correct" match
  mutate(norm_alignment_perc=(alignment_length/max(alignment_length))*100)%>%
  # and filter by the given percentage threshold
  filter(norm_alignment_perc > as.numeric(threshold)) %>%
  # and a resonable e-value
  filter(evalue < .000001) %>%
  # and  count how many good hits remain
  mutate(count = n()) %>%
  # and finally make the dataframe cleaner for plotting by removing duplicates
  filter(!duplicated(count))
  
lineplot <- ggplot(data2, aes(x=flank, y=count, color=region))+
  facet_wrap(~region,ncol = 4, drop = TRUE, scales = "free_x")+
  geom_point()+ scale_x_continuous(breaks=seq(0, 1000, 100))+
  scale_y_continuous(breaks=seq(0, 7, 1))+
  geom_line()+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", size = .1),
    strip.background = element_rect(colour = NA),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.line = element_line(color = 'black', size = .1),
    legend.position=c(0.9, .2)) +
  coord_cartesian(ylim=c(0,7))+
  labs(color="rDNA Region", x="Length of Flanking Region (bp)",
       y="Number of BLAST hits", 
       title="BLASTn Results for BA000007.2 rDNA",
       subtitle=paste0("(Filtered to exclude matches less than ", 
                       threshold,
                       "% of \n query length and hits with E-value >10e-6)"))

# jitter, cause this all the lines converge
data2$countj <- jitter(data2$count, amount = .08)
grouped_lineplot <- ggplot(data2, aes(x=flank, y=countj, color=region))+
  # facet_wrap(~region,ncol = 4, drop = TRUE, scales = "free_x")+
  geom_point()+ 
    scale_x_continuous(breaks=seq(0, 1000, 100), expand = c(0,1))+
  # geom_smooth(se = F)+
      scale_y_continuous(breaks=seq(0, 7, 1))+
  geom_line()+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color=NA, size = .1),
    strip.background = element_rect(colour = NA),
    axis.text.x = element_text(angle=45, hjust=1),
    legend.position=c(0.8, .5),
    axis.line.x= element_line(color = "black", size=.5),
    axis.line.y= element_line(color = "black", size=.5),
    axis.line = element_line(color = NA, size = .1))+
  coord_cartesian(ylim=c(0,7))+
  labs(color="rDNA Region", x="Length of Flanking Region (bp)",
       y="Number of BLAST hits", 
       title="BLASTn Results for BA000007.2 rDNA",
       subtitle=paste0("(Filtered to exclude matches less than ", 
                       threshold,
                       "% of \n query length and hits with E-value >10e-6)"))

pdf(file=file.path(outdir, "sakai_BLAST_results.pdf"), width = 8, height = 6)
lineplot
dev.off()
png(file=file.path(outdir, "sakai_BLAST_results.png"), res=200,  units='in',width = 8, height = 6)
lineplot
dev.off()

pdf(file=file.path(outdir, "grouped_sakai_BLAST_results.pdf"), width = 4, height = 5 )  
grouped_lineplot
dev.off()
png(file=file.path(outdir, "grouped_sakai_BLAST_results.png"), res=200,  units='in',width = 4, height = 5)
grouped_lineplot
dev.off()

