#!/usr/bin/env Rscript

# Install Packages if Needed ----------------------------------------------

# Helper function to install CRAN packages if missing
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Helper for GitHub packages
install_github_if_missing <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools", repos = "https://cloud.r-project.org")
    }
    devtools::install_github(repo)
  }
}

# Install CRAN package
install_if_missing("tidyr")
install_if_missing("reshape2")
install_if_missing("data.table")
install_if_missing("doParallel")
install_if_missing("scales")
install_if_missing("stringr")

# Install GitHub package (replace with actual repo)
install_github_if_missing("cybrBSA", "cbuzby/cybrBSA")

#MAKE FILES FOR ALL CSS EXPERIMENTS FROM OUTPUT TABLE

require(tidyr)
require(cybrBSA)
require(reshape2)
require(data.table)
require(doParallel)
require(scales)
require(stringr)

# Load in Data from command line  ------------------------------------------------------------

myargs <- commandArgs(trailingOnly = TRUE) 

myfile <- myargs[1]
myfilename <- myargs[2]
#myMQCfile <- myargs[3] #don't actually need for this section

#MQC <- read.csv(myMQCfile)

# Load in Data in R ------------------------------------------------------------
#myfile <- "C://Users//cassa//Documents//SiegalLab//Data//Jan25_CSS15.SortedCat.vcf.output.table" #ChangeThis
#MQC <- read.csv("C://Users//cassa//Documents//SiegalLab//Subsampling//MQC_Annotated_REVISED.csv")

#MQC %>% dplyr::select(Pool, PoolID, Dataset = Library, Parent, Bulk, Rep) -> MQC_Key


#GET RID OF EXCLUDED ONES
read.table(myfile, header = TRUE) -> testrawdata
testrawdata %>% dplyr::select(CHROM, POS, ALT) %>% separate(ALT, into = c("First", "Second"), sep = ",") %>% 
  na.omit() -> exclude

ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                       "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                       CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                 "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                 "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                 "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))

exclude %>% dplyr::left_join(.,ChromKey) %>% dplyr::select(-CHROM) %>% dplyr::mutate(CHROM = chromosomes) %>% 
  dplyr::select(-chromosomes) -> exclude

exclude %>% dplyr::mutate(locus = paste(CHROM, POS)) -> exclKEY

cybrInputGATKTable2(myfile) %>% dplyr::mutate(Coverage = as.numeric(AD.REF) + as.numeric(AD.ALT)) %>%
  dplyr::select(POS, CHROM, Dataset, GQ, AD.REF, AD.ALT, Coverage) -> rawdata

parentSNPids <- cybrConvertParentalAlleles(Truncate = TRUE)

rawdata %>% merge(parentSNPids) %>% dplyr::mutate(REFW = as.numeric(Type == "Wine"), REFO = as.numeric(Type == "Oak")) %>%
  group_by(Dataset, CHROM, POS) %>%
  dplyr::mutate(Wine = max(REFW * as.numeric(AD.REF), REFO * as.numeric(AD.ALT)),
         Oak = max(REFW * as.numeric(AD.ALT), REFO * as.numeric(AD.REF))) %>%
  dplyr::select(Dataset, POS, CHROM, Coverage, Wine, Oak) %>%
  pivot_longer(c(Oak, Wine), names_to = "Allele", values_to = "Reads") -> rawdata_called

rawdata_called %>% dplyr::mutate(locus = paste(CHROM, POS)) %>% filter(locus %in% exclKEY$locus == FALSE) %>% 
  dplyr::select(-locus) -> rawdata_called_ex

#saveRDS(rawdata_called_ex, file = "CSS15_called_ex.rds")
finalname_ex <- paste(myfilename, "_called_ex.rds", sep = "")
saveRDS(rawdata_called_ex, file = finalname_ex)

# Smooth Data -------------------------------------------------------------

rawdata_called_ex %>% ungroup() %>% group_by(Dataset) %>% summarize(AvgCovg = mean(Coverage, na.rm = TRUE)) -> Means

rawdata_called_ex %>% dplyr::left_join(Means) %>% ungroup() %>% na.omit() %>%
  group_by(Dataset, CHROM, POS) %>% 
  dplyr::select(Dataset, CHROM, POS, Coverage, AvgCovg) %>%
  distinct() %>%
  dplyr::mutate(CovChange = log(Coverage/AvgCovg)) %>%
  filter(CovChange > -1) -> rawdata_called_ex_Unsmoothed_cf

rawdata_called_ex %>% merge(rawdata_called_ex_Unsmoothed_cf) %>%
  ungroup() %>%
  group_by(CHROM, Allele, Dataset) %>%
  arrange(POS) %>%
  reframe(POS = POS, 
          SmoothCount = frollapply(Reads, n = 200, FUN = cybr_weightedgauss, align = "center")) -> CSS15_G200

#saveRDS(CSS15_G200, file = "CSS15_G200.rds")
finalname <- paste(myfilename, "G200.rds", sep = "")
saveRDS(CSS15_G200, file = finalname)

