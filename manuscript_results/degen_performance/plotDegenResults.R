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
src_folder <-"~/GitHub/riboSeed/manuscript_results/degen_performance/2017-05-16_degenerate_output/genomes/"
files<- dir(src_folder, pattern = ".fasta$")
out_folder <- "~/GitHub/riboSeed/manuscript_results/degen_performance/2017-05-16_results/"
dir.create(out_folder)


for(f in files){
  file.copy(paste0(src_folder, f), paste0(out_folder,f))
  cmd <-paste0("sed -i \"1 s/.*/>", f, "/\" ", out_folder, f)
  print(cmd)
  system(cmd)
}
system(paste0("cat ", out_folder, "0*.fasta > ", out_folder, "combined.fa" ))
mySequences <- readDNAMultipleAlignment(paste0(out_folder, "combined.fa"))
mySequences

degenAlign <- msaConvert(mySequences, type = "seqinr::alignment")
distm <- dist.alignment(degenAlign)
as.matrix(distm)[2:5, , drop=FALSE]
untree<-nj(distm)
plotTree(untree, node.numbers=T)
tree <- root(untree, outgroup = "0.0_reference.fasta", resolve.root = T)
plotTree(tree, node.numbers=T)
plot(tree)
tree2 <- phytools::rotate.multi(tree, 13)
# this needs to be rotated when viewing with plot or plotTree, but not ggtree?
#tree2 <- phytools::rotate.multi(tree2, 21)
plotTree(tree2, node.numbers=T)
tree$tip.label <- gsub("(.*)_(.*)", "\\1_mutation_freq", tree$tip.label)

treedf<-data.frame(name=tree$tip.label,# ord=1:length(tree$tip.label), 
                   # this line gets the actual order its plotted in
                   dist=fortify(tree)[fortify(tree)$isTip==T, "y"])
tree2df<-data.frame(name=tree2$tip.label,# ord=1:length(tree$tip.label), 
                   # this line gets the actual order its plotted in
                   dist=fortify(tree2)[fortify(tree2)$isTip==T, "y"])

(p2<-ggtree(tree2, ladderize = T) + geom_tiplab()
  +geom_text2(aes(subset=!isTip, label=node)))
# 
# (p<-ggtree(tree2, branch.length = "rate", ladderize = T) + 

#   cat ./2017-05-16_degenerate_output/log.txt | grep Changed                         
mutresultraw<-"2017-05-16 16:10:46 - INFO - Changed 850 of 105549 bases
2017-05-16 16:13:32 - INFO - Changed 1342 of 105549 bases
2017-05-16 16:16:05 - INFO - Changed 1930 of 105549 bases
2017-05-16 16:18:29 - INFO - Changed 2327 of 105549 bases
2017-05-16 16:21:00 - INFO - Changed 2654 of 105549 bases
2017-05-16 16:23:32 - INFO - Changed 4251 of 105549 bases
2017-05-16 16:26:04 - INFO - Changed 5947 of 105549 bases
2017-05-16 16:28:49 - INFO - Changed 7329 of 105549 bases
2017-05-16 16:31:09 - INFO - Changed 8441 of 105549 bases
2017-05-16 16:33:41 - INFO - Changed 18987 of 105549 bases
2017-05-16 16:36:29 - INFO - Changed 0 of 105549 bases"

whiraw<-"Processing mutation frequency of 0.0001
Processing mutation frequency of 0.00025
Processing mutation frequency of 0.0005
Processing mutation frequency of 0.00075
Processing mutation frequency of 0.001
Processing mutation frequency of 0.0025
Processing mutation frequency of 0.005
Processing mutation frequency of 0.0075
Processing mutation frequency of 0.01
Processing mutation frequency of 0.05
Processing mutation frequency of 0.0
"
mutresult <-
    gsub("(.* Changed )(\\d*)( .*)", "\\2", strsplit(mutresultraw, "\n")[[1]])

whi <-
  gsub("(.* of )(\\d\\.\\d*)", "\\2", strsplit(whiraw, "\n")[[1]])

dput(tree$tip.label) 
 
nsuc <- read.csv("~/GitHub/riboSeed/manuscript_results/degen_performance/2015-05-aggregated_results.csv", stringsAsFactors = F)
nsuc<-merge(nsuc, treedf, by="name")
(p<-ggtree(tree, ladderize = T) +
    geom_tiplab( align=TRUE, 
                # comment out next line to double check correct labeling
                mapping=aes(label=""), 
                linesize=.5, hjust = -0.05)+
    # ggplot2::xlim(0, 0.5)+
    labs(x="Neighbor-Joining Tree "))

(q <- ggplot(nsuc, aes(y=reorder(name, dist), x=val ))+
    # geom_point(stat="identity", width = .80)+
    geom_vline(0, xintercept = 0, color="grey30") + 
    geom_boxploth(outlier.colour =  "grey50")+
    geom_point(position=position_jitter(width = .05, height = .1))+
    scale_y_discrete(expand=c(.05,.05))+
    scale_x_continuous(expand=c(.1,.1), limits = c(-.1,7), breaks = 0:7)+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(), 
          plot.background = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 20))+
    labs(x="Correctly-assembled rDNAs", title=""))

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
    labs(x="", title=""))


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
