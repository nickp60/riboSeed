# source("http://www.bioconductor.org/biocLite.R")
# biocLite("msa")
# biocLite("ggtree")
library(ggplot2)
library(ape)
library(ggtree)
require(msa)
library(seqinr)
require(phytools)
require(ggstance)
src_folders <-"~/GitHub/riboSeed/manuscript_results/degen_performance"
dir(src_folders)
folders<- dir(src_folders, pattern = "degenerate_output_",full.names = T)
folders

out_folder <- "~/GitHub/riboSeed/manuscript_results/degen_performance/2017-05-17_results/"
dir.create(out_folder)
seeds <- c()
for (d in folders){
(  files <- dir(file.path(d, "genomes"), pattern="*.fasta$", full.names = F))
  files <- gsub("\\.fasta", "", files)
  seed <-gsub("(.*)_(.*?)$", "\\2",  d)
  seeds<-c(seeds, seed)
  for(f in files){
    cmd <-paste0("sed \"1 s/.*/>", f,  "/\" ", d, "/genomes/", f, ".fasta > ", 
                 out_folder, "/", f, "_", seed, ".fasta")
    print(cmd)
    system(cmd)
  }
  system(paste0("cat ", out_folder, "*_", seed, ".fasta > ", out_folder, seed, "_combined.fa" ))
}

treelist <- vector("list", length(seeds))
for (s in 1:length(seeds)){
  print(paste("seed:", seeds[s]))
  mySequences <- readDNAMultipleAlignment(paste0(out_folder,seeds[s], "_combined.fa"))
  print(mySequences)
  degenAlign <-msaConvert(mySequences, type = "seqinr::alignment")

  distm <- dist.alignment(degenAlign)
  untree<-nj(distm)
  tree <- root(untree, outgroup = "0.0_reference", resolve.root = T)
  #plotTree(untree, node.numbers=T)plotTree(tree, node.numbers=T)
  plot(tree)
  treelist[[s]] <- tree
}
class(treelist) <- "multiPhylo"
p<-ggtree(treelist, color="grey70")

tree <- treelist[[2]] 
plotTree(tree)
tree$tip.label <- gsub("(.*)_(.*)", "\\1_mutation_freq", tree$tip.label)

treedf<-data.frame(name=tree$tip.label,# ord=1:length(tree$tip.label), 
                   # this line gets the actual order its plotted in
                   dist=fortify(tree)[fortify(tree)$isTip==T, "y"])

dput(tree$tip.label) 
 
nsuc <- read.csv("~/GitHub/riboSeed/manuscript_results/degen_performance/2015-05-aggregated_results.csv", stringsAsFactors = F)
nsuc<-merge(nsuc, treedf, by="name")
dt<-fortify(tree)
(p<-ggtree(tree, ladderize = T) +
    geom_tiplab( align=TRUE, 
                # comment out next line to double check correct labeling
                mapping=aes(label=""), 
                linesize=.5, hjust = -0.05)+
    # ggplot2::xlim(0, 0.5)+
    labs(x="\n  ", title=""))

(q <- ggplot(nsuc, aes(y=reorder(name, dist), x=val ))+
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
    labs(x="Correctly-assembled rDNAs", title=""))
# multiplot(plotlist = list(p,q), ncol=2, widths = c(.55,.45))

(lab<-ggplot(nsuc, aes(y=reorder(name, dist),
                       x=1,
                       label=gsub("(.*)_mu(.*)","\\1", reorder(name, dist))))+
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
    labs(x="\n", title=""))


multiplot(plotlist = list(p,lab,q), ncol=3, widths = c(.5,.1,.45))

 # d1 <- data.frame(name=tree$tip.label, location=tree$edge.length[1:9])
# p1 <- p %<+% d1 + geom_tippoint(aes(color=location))
# p1
# (p2 <- facet_plot(p, panel="test", data=d1, geom=geom_text, aes(label=d1$name, x=0), color='firebrick'))
# (p3 <- facet_plot(p2, panel="bar", data=nsuc, geom=geom_barh, stat="identity",
#                  mapping=aes(x=val), color='firebrick') + #ggplot2::xlim(0, 7) + theme_tree2()+
#     theme(strip.background = element_blank(), strip.text = element_blank()))
#   
# 
