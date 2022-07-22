#!/usr/bin/env Rscript

start <- Sys.time()

#Load packages
library(ggplot2)
library(tidyr)
library(tidyverse)
library(reshape2)
library(cowplot)
library(foreach)
library(doParallel)

########################################################################################################################
# Load function
########################################################################################################################
pois_binom_sample <- function(lrP, chr = "II", windowsize = 100, coverage = 15634092/2, seedt = Sys.time()){
  lrP <- lrP[lrP$CHROM == chr,]
  
  set.seed(seedt)
  
  #calculate reads by Poisson
  reads <- rpois(n = length(unique(lrP$id)), lambda = mean(lrP$value)/2)
  readsdf <- data.frame(id = unique(lrP$id), numreads = reads)
  lrP %>% left_join(readsdf, by = "id") -> lrP
  
  #calculate binom for number of reads at position
  for(i in 1:length(lrP$CHROM)){
    lrP$subs[i] <- rbinom(1, lrP$numreads[i], prob = lrP$perc[i]) + 1
  }
  #Run the glm on that chromosome
  Results <- foreach (i=unique(lrP$POS), .combine=rbind) %dopar% {
    res <- glm(reads ~ bulk*parent, weights = subs, family = binomial, data = lrP[lrP$POS == i,])
    c(chr, i, res$coefficients)
  }
  
  #Format the results for the next function
  Results <- as.data.frame(Results)
  colnames(Results) <- c("CHROM", "POS", "Intercept", "bulkHigh", "parentWine", "Interaction")
  Results[,2] <- as.numeric(Results[,2])
  Results[,3] <- as.numeric(Results[,3])
  Results[,4] <- as.numeric(Results[,4])
  Results[,5] <- as.numeric(Results[,5])
  Results[,6] <- as.numeric(Results[,6])
  
  Results %>% arrange(POS) %>% select(-Intercept) -> Results
  
  
  #Run the windows on that chromosome
  WResult <- foreach(i=windowsize+1:(length(Results$POS) - (2*windowsize)), .combine=rbind) %dopar% {
    #save the median value from original index 1 to original end index
    res_int <- mean(Results$Interaction[(i-windowsize):(i+windowsize)])
    res_bulk <- mean(Results$bulkHigh[(i-windowsize):(i+windowsize)])
    res_parent <- mean(Results$parentWine[(i-windowsize):(i+windowsize)])
    #print CHROM, index, POS, and median
    c(chr, i, Results$POS[i], Results$bulkHigh[i], Results$parentWine[i],Results$Interaction[i],
      res_bulk, res_parent, res_int)}
  WResult <- as.data.frame(WResult)
  colnames(WResult) <- c("CHROM", "Index", "POS", "bulk_Z", "parent_Z", "Interaction_Z",
                         "bulk_Zprime", "parent_Zprime", "Interaction_Zprime")
  WResult[,2] <- as.numeric(WResult[,2])
  WResult[,3] <- as.numeric(WResult[,3])
  WResult[,4] <- as.numeric(WResult[,4])
  WResult[,5] <- as.numeric(WResult[,5])
  WResult[,6] <- as.numeric(WResult[,6])
  WResult[,7] <- as.numeric(WResult[,7])
  WResult[,8] <- as.numeric(WResult[,8])
  WResult[,9] <- as.numeric(WResult[,9])
  
  return(WResult)
}

########################################################################################################################
# Load data
########################################################################################################################

load("lrPbulksums_p.Rdata")

########################################################################################################################
# Run function to sample
# NEED TO RENAME THESE FOR ACTUAL ARGUMENT - can try on the hpc with this first just to do VII with more samples
########################################################################################################################

pb_VII <- list()
for(i in 1:100){ #PICK NUMBER OF SAMPLES
  pb_VII[[i]] <- pois_binom_sample(lrPbulksums_p, chr = "VII", seedt = i) #PICK CHROMOSOME
}

pb_VII_df <- bind_rows(pb_VII, .id = "Index") #PICK NAME
save(pb_VII_df, file = "pb_VII_df.Rdata") #PICK FILE NAME

print(Sys.time() - start)
