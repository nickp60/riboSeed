
require(dplyr)
require(reshape2)
require(ggplot2)

data <- read.csv2("~/GitHub/riboSeed/manuscript_results/sim_combined.csv", stringsAsFactors = F, sep = ',')

data <- data[data$metric != '# unaligned contigs',]

data$val <- as.numeric(data$val)

# get rid of data we aren't interested in
contig_data <- data[grepl("# contigs \\(", data$metric) & data$cov,]

# I haven't used ggplot in a while, so this is how I choose to rename facets
# I hate it, you hate it -- lets move on
contig_data$bug <- ifelse(grepl("coli", contig_data$bug),"E. coli", "K. pneumoniae")

ggplot(contig_data, 
       aes(x=reorder(metric, as.numeric(gsub("(.*>=)\\s(.*)\\s(.*)", "\\2", contig_data$metric))), y=val, shape=type,color = bug))+
  geom_point()+
  theme(axis.text.x = element_text(angle=320, hjust = 0))+
  facet_wrap(bug~cov, scales="free_y")
fig <- ggplot(contig_data, 
       aes(x=reorder(as.character(cov), cov),
           y=val,
           group=interaction(reorder(metric, as.numeric(gsub("(.*>=)\\s(.*)\\s(.*)", "\\2", contig_data$metric))), type),
          color=reorder(metric, as.numeric(gsub("(.*>=)\\s(.*)\\s(.*)", "\\2", contig_data$metric)))
           ))+
  geom_line(aes(
    linetype=type)) +
  geom_point(aes(shape=type, color=reorder(metric, as.numeric(gsub("(.*>=)\\s(.*)\\s(.*)", "\\2", contig_data$metric)))))+
  scale_y_log10(breaks=c(5,10,20,30,50,100,1000))+
  theme(axis.text.x = element_text(angle=320, hjust = 0))+
  theme_bw()+
  facet_wrap(~bug,scales="free_y")+
  labs(title="Contig Number by Coverage",
       y="Number of Contigs (Log-scaled)",
       x="Simulated Read Fold Coverage",
       color="",
       linetype="Assembly Method")

pdf("./contigs_by_coverage.pdf", width = 6, height = 8)
fig
dev.off()
png(file="./contigs_by_coverage.png",res=200,  units='in',width = 6, height = 8)
fig
dev.off()

