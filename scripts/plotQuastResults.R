
require(dplyr)
require(reshape2)
require(ggplot2)

data <- read.csv2("~/GitHub/riboSeed/manuscript_results/sim_combined.csv", stringsAsFactors = F, sep = ',')

data <- data[data$metric != '# unaligned contigs',]

data$val <- as.numeric(data$val)

contig_data <- data[grepl("# contigs \\(", data$metric),]

ggplot(contig_data, 
       aes(x=reorder(metric, as.numeric(gsub("(.*>=)\\s(.*)\\s(.*)", "\\2", contig_data$metric))), y=val, color = bug))+
  geom_point()+
  theme(axis.text.x = element_text(angle=320, hjust = 0))+
  facet_wrap(~cov, scales="free_y")
