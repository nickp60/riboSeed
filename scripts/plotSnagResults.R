#! /usr/bin/Rscript
#############################################################################################
DEBUG=F
source("~/GitHub/open_utils/NargParse/NargParse.R")
#############################################################################################
help<-"
Blast  Parser
20170103

USAGE:
REQUIRED ARGUMENTS:
[-r resultsfromribosnag.csv ] 
[-o output_directory]	    path to output directory. 
OPTIONAL ARGUMENTS:
[-t <1-100>] bit score percentage similarity to call a hit
[-verb T or F]	          Turn verbose output on and off
"

require(dplyr)
require(reshape2)
require(ggplot2)
this_version<-'0.2'
test_args<-c(
#  "-r", "~/GitHub/riboSeed/sakai_snag_2/BLAST_results/20170103_results_merged.csv",
  "-r", "~/GitHub/riboSeed/manuscript_results/simulated_genome/toyGenome/snag/BLAST_results/20170623_results_merged.csv",
  "-t","75", # threshold percentage
  "-o","~/GitHub/riboSeed/BLAST_results/",  # output directory for figs etc
  "-verb", "T") #Verbose Bool
if (DEBUG | length(commandArgs(T))==0){
  print("Caution! Using with debug dataset")
  pre_args<-test_args
} else{
  pre_args<-commandArgs(T)
}
## @knitr part1

# parse the commandline arguments 
cline_args<-parse_flagged_args(x = pre_args, required = c("-r", "-o"), 
                               optional = c("-t", "-verb", "-p"),
                               version = this_version, help_message = help, test_args=test_args,DEBUG=F)
VERBOSE_bool<-ifelse(arg_bool(cline_args['verb']), T,F)

blast_results<- arg_file(cline_args['r'])




perc_similarity_thresh<-arg_numeric(cline_args["t"], default = 50)*.01
if(perc_similarity_thresh>1 |perc_similarity_thresh<=0){stop("percentage must be between 1 and 100")}
datetime<-gsub("[:-]","",gsub("([\\d:-]*)(\\D*)","\\1",Sys.time()))

output_directory<-arg_directory(cline_args["o"], makedir = T, need_write = T, recursive = F)
output_pre<-paste0(datetime,"_blast_virulence_parser_output_")
setwd(output_directory)
##########
#read in csv
data<-read.csv2(blast_results, header=F, stringsAsFactors = F, sep="\t", 
                col.names = c("query_id", "subject_id", "identity_perc", "alignment_length",
                          "mismatches", "gap_opens", "q_start", "q_end", "s_start",
                          "s_end", "evalue", "bit_score"))
                # colClasses = c("character", "character", "numeric", "integer",
                #                "integer", "integer", "integer", "integer", "integer",
                #                "integer", "numeric", "integer"))
data[,3:ncol(data)] = apply(data[,3:ncol(data)], 2, function(x) as.numeric(as.character(x)))
data$flank <-as.numeric(gsub("(.*)_(\\d*)-bp_flanks", "\\2", data$query_id))
data$region <-(gsub("(.*?)_(.*?)_(.*?)_.*", "\\3", data$query_id))

# make column for ifelse  value/maxvalue>threshold, or if it is threshold or greater than the max bit score per gene.
# data2<-data %>%
#   group_by(query_id) %>%
#   mutate(norm_bit_score=(bit_score/max(bit_score))*100)%>%
#   mutate(pass=ifelse((bit_score/max(bit_score))>perc_similarity_thresh, 1,0))%>%
#   filter(pass==1)

data2<-data %>%
  group_by(query_id) %>%
  mutate(norm_alignment_perc=(alignment_length/max(alignment_length))*100)%>%
  filter(norm_alignment_perc > 90) %>%
  filter(evalue < .000001) %>%
  mutate(count = n()) %>%
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
       title="BLAST Results for BA000007.2 rDNA",
       subtitle="(Filtered to exclude matches less than 90% of \n query length and hits with E-value >10e-6)")

pdf(file="./sakai_BLAST_results.pdf", width = 8, height = 6)
lineplot
dev.off()
png(file="./sakai_BLAST_results.png",res=200,  units='in',width = 8, height = 6)
lineplot
dev.off()
