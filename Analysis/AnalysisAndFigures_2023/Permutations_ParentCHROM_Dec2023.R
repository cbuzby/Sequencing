#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


#Load packages
require(ggplot2)
require(tidyr)
require(dplyr)

#Load Functions
glm_cb2_short <- function(..., W, formula, numgroups = FALSE, outputlength = 4, return = c("Z")) {
  data <- list(...)
  
  #Ensure that there is a formula and W parameter
  if (is.null(W) || is.null(formula)) {
    stop("Weights (W) and formula must be provided")
  }
  #Set formula
  glm_formula <- as.formula(formula)
  #Ensure that formula works for the data provided
  if (!all(names(data) %in% all.vars(glm_formula))) {
    stop("One or more variables in the formula are not provided as arguments")
  }
  
  #########################
  #MAYBEWORKS
  for(i in all.vars(glm_formula)){
    if(length(unique(as.data.frame(data)[,i])) < 2){
      output <- rep(NA, outputlength)
      #print("Not enough levels within groups")
      
      return(output)
    }
  }
  
  glm_fit <- glm(glm_formula, data = as.data.frame(data), weights = W, family = binomial)
  
  #Return specific ones like z-score only
  if(return %in% "Z"){
    output <- summary(glm_fit)$coefficients[((length(summary(glm_fit)$coefficients)*0.5)+1):((length(summary(glm_fit)$coefficients)*0.75))]
  }
  
  if(length(output) == outputlength){
    return(output)
  }else{
    return(rep(NA, outputlength))
  }
  
}


# Replicates, fix CHROM, Parent
cybrPermute_byCHRParent_Rep <- function(dataset){
  #start.time <- Sys.time()
  
  #print("Make sure that dilute bulk is labeled D")
  
  dataset %>% 
    distinct() %>% ungroup() %>%
    group_by(CHROM, POS, Allele, 
             Bulk, 
             Rep, #might not have these
             Parent) %>% 
    summarize(culprits = length((SmoothCount))) %>% 
    merge(dataset) %>% 
    filter(culprits == 1) %>% 
    ungroup() %>%
    distinct() %>% #THIS IS IMPORTANT
    pivot_wider(names_from = Allele, values_from = SmoothCount) -> newnewtest
  
  #these are now all of the ones that can be used to permute
  newnewtest %>% filter(Bulk == "D") %>% select(CHROM, POS, Parent, Oak, Wine) %>% 
    group_by(CHROM, Parent) %>%
    summarize(CHROM = CHROM, Parent = Parent, Oak = Oak, Wine = Wine, POS = POS,
              POS2 = sample(POS)) %>%
    group_by(CHROM, Parent, POS2) %>%
    summarize(CHROM = CHROM, Parent = Parent, Oak = Oak, Wine = Wine,
              POS = POS2, 
              Rep = sample(c("a", "b"))) %>%
    mutate(Bulk = "A") -> shuffled_DiluteA
  
  newnewtest %>% filter(Bulk == "D") %>% select(CHROM, POS, Parent, Oak, Wine) %>% 
    group_by(CHROM, Parent) %>%
    summarize(CHROM = CHROM, Parent = Parent, Oak = Oak, Wine = Wine, POS = POS,
              POS2 = sample(POS)) %>%
    group_by(CHROM, Parent, POS2) %>%
    summarize(CHROM = CHROM, Parent = Parent, Oak = Oak, Wine = Wine,
              POS = POS2, 
              Rep = sample(c("a", "b"))) %>%
    mutate(Bulk = "B") -> shuffled_DiluteB
  
  rbind(shuffled_DiluteA, shuffled_DiluteB) %>% pivot_longer(c(Oak, Wine), names_to = "Allele", values_to = "SmoothCount") -> shuffletoglm
  
  #Trying this again
  
  shuffletoglm %>% na.omit() %>%
    #Original Script
    group_by(CHROM, POS) %>%
    mutate_if(is.character, as.factor) %>%
    summarize(Summary = glm_cb2_short(Allele = Allele,
                                      Bulk = Bulk,
                                      Parent = Parent,
                                      Rep = Rep,
                                      W = SmoothCount,
                                      formula = "Allele ~ Bulk * Parent + Rep",
                                      numgroups = 16, outputlength = 5),
              Factor = (c("Intercept", "Bulk", "Parent", "Rep", "Interaction"))) -> glmresult
  #end.time = Sys.time()
  #print(end.time - start.time)
  return(glmresult)
}

################################################################################

mydata <- readRDS(args[1])

cybrPermute_byCHRParent_Rep(mydata) %>% ungroup() %>% group_by(CHROM, Factor) %>% summarize(quant095 = quantile(abs(Summary), 0.95, na.rm = TRUE),
                                                                                  quant099 = quantile(abs(Summary), 0.99, na.rm = TRUE),
                                                                                  quant0995 = quantile(abs(Summary), 0.995, na.rm = TRUE)) %>%
  mutate(array = args[2])-> result

print(as.data.frame(result))
