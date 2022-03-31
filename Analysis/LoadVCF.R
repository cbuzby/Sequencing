## Visualizing the snp.vcf from sequencing pipleine
## 1/17/2022

choose.files()
tmp_vcf<-readLines("C:\\Users\\cassa\\OneDrive\\Documents\\R_SiegalLab\\SequencingPractice\\NZ_UnselectedA.snps.vcf")

NZ_UA_vcg <- tmp_vcf[50:length(tmp_vcf)]
NZ_UA_vcg_unlisted <-unlist(strsplit(NZ_UA_vcg[1:length(NZ_UA_vcg)],"\t"))

NZ_UA_mat <- matrix(NZ_UA_vcg_unlisted, ncol = 10, byrow = TRUE)
NZUA_df <- as.data.frame(NZ_UA_mat)
colnames(NZUA_df) <- NZUA_df[1,]
NZUA_df <- NZUA_df[2:length(NZUA_df$POS),]

head(NZUA_df)
colnames(NZUA_df) <- c("CHROM", "POS", "ID", "REF","ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NZ_UnselectedA.sort")

library(ggplot2)
ggplot(NZUA_df, aes(x = as.numeric(QUAL), color = CHROM)) + geom_density() + ggtitle("Quality by Chromosome") + theme_classic()

##################

# Reformatting to make the QTLseqr input table

head(unique(NZUA_df$REF), 100)

NZUA_df[NZUA_df$REF == "attgttgttgttgttgtt",]
