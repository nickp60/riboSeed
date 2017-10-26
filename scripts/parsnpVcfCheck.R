### Comparing the de fere novo assembly to the true assembly

vcf_path <- "~/GitHub/riboSeed/manuscript_results/pao_results/2017-10-25-pao_ATCC_ref/parsnp_results/parsnp.vcf"
vcf <- read.csv2(vcf_path, skip = 5, sep="\t", stringsAsFactors = F)
head(vcf)

gff_names <- c("seq", "program", "feature", "start", "end", "score", "strand", "unk", "anno")
gff_path <- "~/GitHub/riboSeed/manuscript_results/pao_results/scan/scannedScaffolds.gff"
gff <-read.csv2(gff_path, stringsAsFactors = F, header=F, skip = 1, sep="\t", col.names = gff_names)
str(vcf)

### get SNPs in de fere novo rDNAs
rdnas_coords <- data.frame(start=c(399001, 1040539, 1863045, 2810154) - 1000,
                           end=c(404148, 1045687, 1868194, 2815303) + 1000)
mark_rdna <- function(vcf, col, rdnas_coords){
  # vcf$rdna <- 0
  # for (j in 1:nrow(rdnas_coords)){  # for each rDNA coord set
  #   for (i in 1:nrow(vcf)){  # for each line in vcf
  #     # if POS is within the coord set, change rdna from 0 to index of coord set 
  #     if (vcf[i, "POS"] > rdnas_coords[j, "start"] & vcf[i, "POS"] < rdnas_coords[j, "end"]){
  #       vcf[i, "rdna"] <- paste0("(", rdnas_coords[j, "start"], " - ", rdnas_coords[j, "end"], ")")
  #     }
  #   }
  # }
  vcf <- vcf[col != 0, ]
  vcf$rdna <- 0
  for (j in 1:4){
    cluster <- paste0("(", rdnas_coords[j, "start"], " - ", rdnas_coords[j, "end"], ")")
    vcf$rdna <- ifelse(
      vcf$rdna == 0,
      ifelse( vcf$POS > rdnas_coords[j, "start"] & 
                vcf$POS < rdnas_coords[j, "end"] &
                vcf[, col] != 0, cluster, 0),
      vcf$rdna       
    )
    print(j)
  }
  return(vcf)
}
print("Total SNPs in the de fere novo assembly:")
print(sum(vcf$NZ_CP017149.1_de_fere_novo_contigs.fasta))
print("Total 'SNPs' in the reference:")
print(sum(vcf$NZ_CP017149.fasta))
de_fere_vcf <- mark_rdna(vcf=vcf, col="NZ_CP017149.1_de_fere_novo_contigs.fasta", rdnas_coords = rdnas_coords)
ref_vcf <- mark_rdna(vcf=vcf, col="NZ_CP017149.fasta", rdnas_coords = rdnas_coords)
pretty_table_defere <- de_fere_vcf[de_fere_vcf$rdna != 0, c("rdna", "POS", "REF", "ALT")]
print(pretty_table_defere)
print("SNPs in rDNA regions of de fere novo assembly:")
sum(table(de_fere_vcf$rdna)[2:length(table(de_fere_vcf$rdna))])
print("SNPs in rDNA regions of the reference:")
sum(table(ref_vcf$rdna)[2:length(table(ref_vcf$rdna))])



