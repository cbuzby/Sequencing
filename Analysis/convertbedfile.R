#Convert BED file to the new chromosome names
setwd("C:/Users/cassa/OneDrive/Documents/R_SiegalLab/SequencingPractice/")

ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", 
                                       "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"), 
                       CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10", 
                                 "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6", 
                                 "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5", 
                                 "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))

bed <- read.table("C:/Users/cassa/OneDrive/Documents/R_SiegalLab/SequencingPractice/excluded_regions.bed")


bed %>% mutate(chromosomes = V1) %>% merge(.,ChromKey) %>% transmute(V1 = CHROM, V2 = V2, V3 = V3) -> bed2

write.table(bed2, file = "excludedregions_NC.bed", row.names = FALSE, col.names = FALSE, quote = FALSE)
