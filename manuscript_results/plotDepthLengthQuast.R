# 20170518 NRW
require(ggplot2)
folder <- "~/GitHub/riboSeed/manuscript_results/simulated_reads/"
folder
files <- dir(folder, pattern = "*_cov_*", full.names = T)
files

f <- files[1]
df <- data.frame("metric"=NA, "de_fere_novo"=NA, "de_novo"=NA, "cov"=NA, "len"=NA, "bug"=NA)
cnames <- c("metric","de_fere_novo","de_novo","cov","len", "bug")
for(f in files){
  base <- basename(f)
  cov <- as.numeric(gsub("(.*)_cov_(.*).csv", "\\2", base))
  len <- as.numeric(gsub("(.*)_len_(.*)_cov(.*).csv", "\\2", base))
  bug <- gsub("(.*)_sim_(.*)_cov(.*).csv", "\\1", base)
  data <- read.csv2(f, stringsAsFactors = F, sep = "\t", header = F,col.names = cnames )
  data$cov <- cov 
  data$len <- len
  data$bug <- bug
  df <- rbind(df, data)
    
}

df2 <- df[df$metric %in% c(
  "# contigs (>= 0 bp)"
  # "# contigs (>= 1000 bp)",
  # "# contigs (>= 5000 bp)",
  # "# contigs (>= 10000 bp)",    
  # "# contigs (>= 25000 bp)",
  # "# contigs (>= 50000 bp)"
  ), ]
df2$dif <-as.numeric(df2$de_fere_novo) - as.numeric(df2$de_novo)

df2$metric <- reorder(df2$metric, -as.numeric(gsub("(.*)= (\\d*) bp.*", "\\2", df2$metric)))
df2$metric_num <- as.numeric(gsub("(.*)= (\\d*) bp.*", "\\2", df2$metric))

ggplot(df2, aes(x=reorder(as.character(cov), cov), y=dif, 
                color=as.character(len), 
                fill=as.character(len), 
                shape=as.character(len))) +
  geom_bar(position="dodge", stat="identity", width=.4) +
  facet_wrap(~bug, ncol = 2, scales = "free") +
  # geom_line(aes(group=len)) +
  scale_x_discrete(expand=c(.01,.01)) + #, limits = c(-.05,7.5), breaks = 0:7)+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  # scale_color_discrete(guide=F)+
  labs(x="Coverage Depth",
       y="Contig Number Difference \n(De fere novo - De novo)",
       fill="Read Length",
       color="Read Length",
       shape="Coverage")+
  geom_hline(0, yintercept = 0, color="grey30") + 
  scale_y_continuous(expand=c(.01,.01))+
  theme(#axis.title.y = element_blank(),
        # axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size=15))+
  coord_flip() 

