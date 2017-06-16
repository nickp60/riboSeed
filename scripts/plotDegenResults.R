#!/usr/bin/env Rscript
# version 0.0.2
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("msa")
# biocLite("ggtree")
# require(ggplot2)

require(argparse, quietly = T, warn.conflicts = F)
require(reshape2, quietly = T, warn.conflicts = F)
require(msa, quietly = T, warn.conflicts = F)
require(seqinr, quietly = T, warn.conflicts = F)
require(ape, quietly = T, warn.conflicts = F)
require(ggtree, quietly = T, warn.conflicts = F)
require(phytools, quietly = T, warn.conflicts = F)
require(ggstance, quietly = T, warn.conflicts = F)

###################################################################
# create parser object
parser <- ArgumentParser()
# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-i", "--input_dir", action="store",
                    dest="src_folder", help="root of run_degen_sim.sh output")
parser$add_argument("-r", "--results_csv", action="store", 
                    default="NA",
                    dest="results", help="columns must be name, seed, value, variable")
parser$add_argument("-o", "--out_dir", action="store",
                    dest="out_folder", help="output_directory")

args <- parser$parse_args()
print(args)
if (dir.exists(args$out_folder)){
  print("output dir exists!")
  quit()
  
}

(analysis_date <- gsub("(.*)([0-9]{4}-[0-9]{2}-[0-9]{2})(.*)", "\\2", args$src_folder))

(names <- c("sourcepath", "assembly", "total", "good",  "ambiguous", "bad"))
(dir(args$src_folder))
(reports <- dir(paste0(file.path(args$src_folder, analysis_date), "_reports/")))
print(length(reports))
if(length(reports)==0){
  print("no reports found")
  quit()
}
# (out_folder <- paste0("~/GitHub/riboSeed/", analysis_date, "_degen_results/"))

results_table <- as.data.frame(setNames(replicate(6,numeric(0), simplify = F), names))
# results_table_old <- as.data.frame(setNames(replicate(6,numeric(0), simplify = F), names))
# for (report in reports){
#   try(
#     results_table_old <- rbind(
#       results_table_old, read.csv2(
#         paste0(args$src_folder, analysis_date, "_reports/", report), 
#         sep="\t", stringsAsFactors = F, header = F, col.names = names))
#   )
# }

dir.create(args$out_folder)
dir(args$src_folder)


folders<- dir(args$src_folder, pattern = "_degenerate_output_",full.names = T)
folders

seeds <- c()
ngenome_list <- c()
for (d in folders){
  (files <- dir(file.path(d, "genomes"), pattern="*.fasta$", full.names = F))
  ngenome_list <- c(ngenome_list, files)
  files <- gsub("\\.fasta", "", files)
  seed <-gsub("(.*)_(.*?)$", "\\2",  d)
  seeds<-c(seeds, seed)
  for(f in files){
    cmd <-paste0("sed \"1 s/.*/>", f,  "/\" ", d, "/genomes/", f, ".fasta > ", 
                 args$out_folder, "/", f, "_", seed, ".fasta")
    print(cmd)
    system(cmd)
  }
  system(paste0("cat ", args$out_folder, "*_", seed, ".fasta > ", 
                args$out_folder, seed, "_combined.fa" ))
}
print("Finished rewriting fastas with better accessions and combining into MSAs")


###################################################################
(freqs <- gsub("(.*)_(.*)", "\\1", unique(ngenome_list)))
sources_list <- c()
for (report in reports){
  report_path <- paste0(
    file.path(args$src_folder, analysis_date), "_reports/", report)
  if (file.exists(report_path)){
    new_table <- read.csv2(
      report_path, sep="\t", stringsAsFactors = F, header = F, col.names = names)
#    sources_list <- c(sources_list, unique())
    new_table$freq <- gsub("(.*)_(.*?)\\.fasta", "\\2", new_table$assembly)
    head_name <- gsub("(.*)_(.*?)\\.fasta", "\\1", new_table$assembly)
    print(nrow(new_table))
    for (freq in freqs){
      if (freq %in% new_table$freq){
        print("found")
      } else{
        print(paste("missing", freq))
        new_table <- rbind(
          new_table, 
          data.frame("sourcepath" = new_table$sourcepath[1], 
                     "assembly" = paste0(head_name[1], "_", freq, ".fasta"), 
                     "total" = max(new_table$total),
                     "good" = 0,
                     "ambiguous" = 0,
                     "bad" = 0,
                     "freq" = freq)
          )        
      }
    }
    new_table$freq <- NULL
    results_table <- rbind(results_table, new_table)
    print("\n\n")
    print(results_table)
    print("\n\n")
    
  }
}
str(results_table)
results_table$freq <- NULL
results_table$name <-as.numeric(
  gsub("(.*)_(.*)", "\\2", gsub("\\.fasta", "", results_table$assembly)))
str(results_table)
results_table$assembly <- NULL
tall <- melt(results_table, 
             id.vars = c("sourcepath", "name", "total"), 
             measure.vars = c("good", "bad", "ambiguous"), factorsAsStrings = T)
str(tall)
###################################################################
if (args$results != "NA"){
  oldtall <- tall
  print("Using Pre-computed scoring results")
  tall <- read.csv(args$results, stringsAsFactors = F)
  tall$sourcepath <- oldtall[oldtall$variable=="good", "sourcepath"]
  tall$seed <-NULL
  tall$variable <- as.factor(tall$variable)
  str(tall)
} else {
  print("no")
}

treelist <- vector("list", length(seeds))
print("Processing each MSA")
for (s in 1:length(seeds)){
  print(paste("seed:", seeds[s]))
  mySequences <- readDNAMultipleAlignment(paste0(args$out_folder,seeds[s], "_combined.fa"))
  # print(mySequences)
  degenAlign <-msaConvert(mySequences, type = "seqinr::alignment")

  distm <- dist.alignment(degenAlign)
  untree<-nj(distm)
  tree <- root(untree, outgroup = "0.0_reference", resolve.root = T)
  #plotTree(untree, node.numbers=T)plotTree(tree, node.numbers=T)
  #plot(tree)
  treelist[[s]] <- tree
}
class(treelist) <- "multiPhylo"
p<-ggtree(treelist, color="grey70")

tree <- treelist[[2]] 
# plotTree(tree)

tree$tip.label <- gsub("(.*)_(.*)", "\\1_mutation_freq", tree$tip.label)

treedf<-data.frame(label=tree$tip.label,
                   name=as.numeric(gsub("(.*?)_.*", "\\1", tree$tip.label)),
                   # this line gets the actual order its plotted in
                   dist=fortify(tree)[fortify(tree)$isTip==T, "y"])

#dput(tree$tip.label) 

print(str(treedf ))
print(str(tall))
print("merging riboScore results with tree DF" )
nsuc<-merge(tall, treedf, by="name")
print(nsuc)
write.csv(nsuc, file = file.path(args$out_folder, "results.csv"))
dt<-fortify(tree)
p<-ggtree(tree, ladderize = T) +
    geom_tiplab( align=TRUE, 
                # comment out next line to double check correct labeling
                mapping=aes(label=""), 
                linesize=.5, hjust = -0.05)+
    # ggplot2::xlim(0, 0.5)+
    labs(x="\n  ", title="")

q <- ggplot(nsuc[nsuc$variable == "good",],
             aes(y=reorder(name, dist), x=value ))+
    # geom_point(stat="identity", width = .80)+
    geom_vline(0, xintercept = 0, color="grey30") + 
    geom_boxploth(outlier.colour =  NA,color="grey40", fill="grey60",width=.7)+
    geom_point(position=position_jitter(width = .02, height = .12), shape=1)+
    scale_y_discrete(expand=c(.01,.01))+
    scale_x_continuous(expand=c(.01,.01), limits = c(-.05,7.5), breaks = 0:7)+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(), 
          plot.background = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size =15))+
    labs(x="Correctly-assembled rDNAs", title="")
# multiplot(plotlist = list(p,q), ncol=2, widths = c(.55,.45))
(unique_idx <- nsuc[!duplicated(nsuc$name),])

options(scipen = 99999)
nsuc$label <- as.character(nsuc$name)
print(unique(nsuc$label))
lab<-ggplot(nsuc[!duplicated(nsuc$name),],
             aes(y=reorder(name, dist),
                 x=1,
                 # label=gsub("(.*)_mu(.*)","\\1", reorder(name, dist))))+
                 label=label))+
  theme(
      line = element_blank(), 
      rect = element_blank(), 
      text = element_text(face = "plain", colour = "black", 
                          lineheight = 0.9, hjust = 0.5, 
                          vjust = 0.5, angle = 0, margin = margin()),
      axis.text = element_blank(), 
      axis.title.y = element_blank(), 
      legend.text = element_text(size = rel(0.8)), 
      legend.title = element_text(hjust = 0), 
      strip.text = element_text(size = rel(0.8)) 
      # plot.margin = unit(c(0, 0, 0, 0), "lines")
      )+
    geom_text(size=3, color="grey20")+
    labs(x="\n", title="")


pdf(file = paste0(args$out_folder, "degenerate_flanking_regions_plot.pdf"), width = 8, height = 4)
multiplot(plotlist = list(p,lab,q), ncol=3, widths = c(.5,.1,.45))
dev.off()
png(filename = paste0(args$out_folder, "degenerate_flanking_regions_plot.png"), width = 8, height = 4, units = "in", res = 300)
multiplot(plotlist = list(p,lab,q), ncol=3, widths = c(.5,.1,.45))
dev.off()

quit()
 
