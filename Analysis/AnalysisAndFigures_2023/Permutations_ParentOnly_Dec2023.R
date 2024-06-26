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


# Replicates, fix Parent
cybrPermute_byParent_Rep <- function(dataset){
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
  newnewtest %>% filter(Bulk == "D",
                        CHROM %in% c("I", "III", "V", "VIII", "M") == FALSE) %>% 
    select(CHROM, POS, Parent, Oak, Wine) %>% 
    group_by(Parent) %>% mutate(Loc = paste(CHROM, POS)) %>%
    summarize(Parent = Parent, Oak = Oak, Wine = Wine, 
              Loc = sample(Loc)) %>%
    group_by(Parent, Loc) %>%
    summarize(Parent = Parent, Oak = Oak, Wine = Wine,
              Loc = Loc, 
              Rep = sample(c("a", "b"))) %>%
    mutate(Bulk = "A") -> shuffled_DiluteA2
  
  newnewtest %>% filter(Bulk == "D",
                        CHROM %in% c("I", "III", "V", "VIII", "M") == FALSE) %>% 
    select(CHROM, POS, Parent, Oak, Wine) %>% 
    group_by(Parent) %>% mutate(Loc = paste(CHROM, POS)) %>%
    summarize(Parent = Parent, Oak = Oak, Wine = Wine, 
              Loc = sample(Loc)) %>%
    group_by(Parent, Loc) %>%
    summarize(Parent = Parent, Oak = Oak, Wine = Wine,
              Loc = Loc, 
              Rep = sample(c("a", "b"))) %>%
    mutate(Bulk = "B") -> shuffled_DiluteB2
  
  rbind(shuffled_DiluteA2, shuffled_DiluteB2) %>% pivot_longer(c(Oak, Wine), names_to = "Allele", values_to = "SmoothCount") -> shuffletoglm2
  
  #Trying this again
  
  shuffletoglm2 %>% na.omit() %>%
    #Original Script
    group_by(Loc) %>%
    mutate_if(is.character, as.factor) %>%
    summarize(Summary = glm_cb2_short(Allele = Allele,
                                      Bulk = Bulk,
                                      Parent = Parent,
                                      Rep = Rep,
                                      W = SmoothCount,
                                      formula = "Allele ~ Bulk * Parent + Rep",
                                      numgroups = 16, outputlength = 5),
              Factor = (c("Intercept", "Bulk", "Parent", "Rep", "Interaction"))) -> glmresult2
  
  #end.time = Sys.time()
  #print(end.time - start.time)
  
  return(glmresult2)
}


################################################################################

mydata <- readRDS(args[1])

cybrPermute_byParent_Rep(mydata) %>% ungroup() %>% group_by(Factor) %>% summarize(quant095 = quantile(abs(Summary), 0.95, na.rm = TRUE),
                                                                                  quant099 = quantile(abs(Summary), 0.99, na.rm = TRUE),
                                                                                  quant0995 = quantile(abs(Summary), 0.995, na.rm = TRUE)) %>%
  mutate(array = args[2])-> result

print(as.data.frame(result))

print(vector(result))
