
# BiocInstaller::biocLite("seqinr")
require(ggplot2)
require(seqinr)
# ref <- seqinr::read.fasta("~/GitHub/riboSeed/manuscript_results/simulated_genome/BA000007.2.fasta")
ref_path <- "~/GitHub/riboSeed/manuscript_results/entropy/sakai_1st_rDNA_scan/contigs/BA000007.2_0_226105..233199_0.fa"
ref <- seqinr::read.fasta(ref_path)
# vcf <- read.csv2("~/testparsnp/P_2017_05_31_082508757460/parsnp", skip = 5, sep="\t", stringsAsFactors = F)
vcf_path <- "~/GitHub/riboSeed/2017-10-13-parsnp/parsnp.vcf"
vcf <- read.csv2(vcf_path, skip = 5, sep="\t", stringsAsFactors = F)
head(vcf)
# these are the columns that actually have the snps for each of our samples
(cols <- colnames(vcf)[10:ncol(vcf)])
pre = c()
within = c()
post = c()
throughout = c()
# the length of the reference sequence
ref_len <- length(ref$BA000007.2_0_226105..233199)

# for each sample,
for( i in colnames(vcf)[grepl(".fasta", colnames(vcf))]){
  # we assume that we are working with 1kb flanking
  pre = c(pre, sum(vcf[vcf$POS < 1000, i]) / 1000 )
  #  the within region is the regions without the 1kb on either side
  within  <- c(within, sum(vcf[vcf$POS > 1000 & vcf$POS < (ref_len - 1000) , i]) / (ref_len-2000) )
  # the whole sequence
  throughout  <- c(throughout, sum(vcf[, i]) / (ref_len) )
  # the last 1kb in the sequence 
  post <- c(post, sum(vcf[vcf$POS > (ref_len-1000), i]) / 1000 )
}
mean(pre)
mean(within)
mean(post)
mean(throughout)
write(x = c(paste0("reference:\t", ref_path), 
            paste0("Parsnp vcf:\t", vcf_path),
            "#--- Average Substitution Rates ---",
            paste0("5' 1kb flanking region:\t", mean(pre)),
            paste0("16S, 5S, 23S, and intergenic\t:", mean(within)),
            paste0("Entire Region:\t", mean(throughout))),
      file = "~/GitHub/riboSeed/2017-10-13-parsnp/substitution_results.txt")

plotdf <- data.frame(position=1:ref_len)
vcf$total <- rowSums(vcf[,grepl(".fasta", colnames(vcf))])
vcf$total


plotdf <-merge(plotdf, vcf, by.x="position", by.y="POS", all.x = T)  
plotdf$total[is.na(plotdf$total)] <- 0
ggplot(plotdf, aes(x=position, y=total))+
  geom_point()


