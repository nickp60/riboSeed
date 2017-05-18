# 20170518 NRW
require(ggplot2)
folder <- "~/GitHub/riboSeed/manuscript_results/simulated_reads_coli/"
folder
files <- dir(folder, pattern = "*_cov_*", full.names = T)
files

f <- files[1]
df <- data.frame("metric"=NA, "de_fere_novo"=NA, "de_novo"=NA, "cov"=NA, "len"=NA)
cnames <- c("metric","de_fere_novo","de_novo","cov","len")
for(f in files){
  base <- basename(f)
  cov <- as.numeric(gsub("(.*)_cov_(.*).csv", "\\2", f))
  len <- as.numeric(gsub("(.*)_len_(.*)_cov(.*).csv", "\\2", f))
  data <- read.csv2(f, stringsAsFactors = F, sep = "\t", header = F,col.names = cnames )
  data$cov <- cov 
  data$len <- len
  df <- rbind(df, data)
    
}
df2 <- df[df$metric %in% c(
  "# contigs (>= 0 bp)",
  "# contigs (>= 1000 bp)",
  "# contigs (>= 10000 bp)",  
  "# contigs (>= 25000 bp)",
  "# contigs (>= 5000 bp)",
  "# contigs (>= 50000 bp)" ), ]
df2$dif <-as.numeric(df2$de_fere_novo) - as.numeric(df2$de_novo)

ggplot(df2, aes(x=reorder(metric, -as.numeric(gsub("(.*)= (\\d*) bp.*", "\\2", df2$metric))), y=dif, color=as.character(len))) +
  geom_point() + facet_wrap(~cov, ncol = 4, scales = "free_x") +
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  coord_flip() + 
  labs(x="Contigs", y="De fere novo - De novo", color="Read Length")



