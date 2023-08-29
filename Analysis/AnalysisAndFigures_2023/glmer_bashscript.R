#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#INPUTS: 
#(1) File (.rds) of rolled data with structure CHROM | POS | BulkName_ParentName_RepName_AlleleName | Nextbulk... etc
#(2) Chromosome to filter by, as character

#OUTPUTS:
#(1) File saved as arg1_arg2_glmoutput.rds, which will have the GLMs for each chromosome saved; 
# -> Uses the formula "Allele ~ Bulk*Parent+(1 | Rep)"
# -> Includes ALL coefficients: Estimate, SE, Zscore, P-value

################################################################################

#Load packages
library(tidyr)
library(reshape2)
library(data.table)
library(dplyr)
library(lme4)

#Write the function ############################################################

mixedglm <- function(..., W, formula) {
  data <- list(...)
  
  if (is.null(W) || is.null(formula)) {
    stop("Weights (W) and formula must be provided")
  }
  
  glm_formula <- as.formula(formula)
  
  if (!all(names(data) %in% all.vars(glm_formula))) {
    stop("One or more variables in the formula are not provided as arguments")
  }
  
  #CHANGE THIS TO ADJUST TYPE OF REGRESSION
  #glm_fit <- lme4(glm_formula, data = as.data.frame(data), weights = W, family = binomial)
  
  b <- glmer(glm_formula, 
             weights = W, family = binomial, 
             data = as.data.frame(data), nAGQ=20)
  
  return(as.vector(summary(b)$coefficients))
  
}

#Import Data ###################################################################
rollData <- readRDS(args[1])


#Set the filename ##############################################################
arg1file <- gsub(pattern = ".rds", replacement = "", x = args[1])
glmfilename <- paste(arg1file, args[2], "glmoutput.rds", sep = "_")

#Run the code for the glm ######################################################
rollData %>% filter(CHROM == args[2]) %>%
  pivot_longer(c(-CHROM, -POS), names_to = c("Bulk", "Parent", "Rep", "Allele"), 
                          names_sep = "_", values_to = "Reads") %>% 
  mutate_if(is.character, as.factor) %>%
  group_by(CHROM, POS) %>%
  reframe(GLM = mixedglm(W = Reads, formula = "Allele ~ Bulk*Parent+(1 | Rep)",
                         Allele = Allele, Bulk = Bulk, Parent = Parent, Rep = Rep),
          Effect = rep(c("Intercept","Bulk", "Parent", "Interaction"),4),
          Type = c(rep("Effect",4), rep("SE",4), rep("Zscore", 4), rep("pval", 4))) -> glmoutput


#Save the file #################################################################
saveRDS(glmoutput, file = glmfilename)
