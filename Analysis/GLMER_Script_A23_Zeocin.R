#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Load packages
library(tidyr)
library(reshape2)
library(dplyr)
library(lme4)

#Load in function
glm_mixed_model_Zeo <- function(Dilute_Oak_A_Oak, Dilute_Oak_A_Wine, Dilute_Oak_B_Oak, Dilute_Oak_B_Wine, 
                                Zeocin_Oak_A_Oak, Zeocin_Oak_A_Wine,Zeocin_Oak_B_Oak, Zeocin_Oak_B_Wine,
                                Dilute_Wine_A_Oak ,Dilute_Wine_A_Wine, Dilute_Wine_B_Oak, Dilute_Wine_B_Wine,
                                Zeocin_Wine_A_Oak,Zeocin_Wine_A_Wine, Zeocin_Wine_B_Oak,Zeocin_Wine_B_Wine){
  
  combineddata <- data.frame(Bulk = factor(c("L", "L","L", "L", 
                                             "H", "H","H", "H", 
                                             "L", "L","L", "L", 
                                             "H", "H","H", "H")),
                             Parent = factor(c("O", "O", "O", "O", 
                                               "O", "O", "O", "O",
                                               "W", "W", "W", "W", 
                                               "W", "W", "W", "W")),
                             Allele = factor(c("O", "W", "O", "W",
                                               "O", "W", "O", "W", 
                                               "O", "W", "O", "W",
                                               "O", "W", "O", "W")),
                             Rep = factor(1:16),
                             Reads = c(Dilute_Oak_A_Oak, Dilute_Oak_A_Wine, Dilute_Oak_B_Oak, Dilute_Oak_B_Wine, 
                                       Zeocin_Oak_A_Oak, Zeocin_Oak_A_Wine,Zeocin_Oak_B_Oak, Zeocin_Oak_B_Wine,
                                       Dilute_Wine_A_Oak ,Dilute_Wine_A_Wine, Dilute_Wine_B_Oak, Dilute_Wine_B_Wine,
                                       Zeocin_Wine_A_Oak,Zeocin_Wine_A_Wine, Zeocin_Wine_B_Oak,Zeocin_Wine_B_Wine))
  
  combineddata %>% group_by(Bulk, Parent) %>% summarize(key = factor(paste(Bulk, Parent, sep = "_"))) %>% merge(combineddata) -> combineddata
  
  b <- glmer(Allele ~ Bulk*Parent+(1+key | Rep), 
             weights = Reads, family = binomial, 
             data = combineddata)
  
  
  return(c(summary(b)$coefficients[1:12]))
  
} 


#Load data
Zeo <- readRDS("Zeocin_cybr2.rds")

#Roll data
Zeo %>% pivot_longer(c(-CHROM, -POS), names_to = "label") %>% 
  group_by(CHROM, label) %>% arrange(POS) %>% 
  summarize(POS = POS, 
            CHROM = CHROM, 
            SmoothCount = ceiling(frollapply(value, n = 200, FUN = median))) %>% #CHANGE THIS LINE
  na.omit() %>% 
  pivot_wider(names_from = label,values_from = SmoothCount) -> rollData


#Run function
rollData %>% group_by(CHROM, POS) %>% #head(2) %>%
  summarise(summary = glm_mixed_model_Zeo(Dilute_Oak_A_Oak, Dilute_Oak_A_Wine, Dilute_Oak_B_Oak, Dilute_Oak_B_Wine, 
                                          Zeocin_Oak_A_Oak, Zeocin_Oak_A_Wine,Zeocin_Oak_B_Oak, Zeocin_Oak_B_Wine,
                                          Dilute_Wine_A_Oak ,Dilute_Wine_A_Wine, Dilute_Wine_B_Oak, Dilute_Wine_B_Wine,
                                          Zeocin_Wine_A_Oak,Zeocin_Wine_A_Wine, Zeocin_Wine_B_Oak,Zeocin_Wine_B_Wine),
            label = c("E_intercept", "E_Bulk", "E_Parent", "E_Interaction",
                      "SE_intercept", "SE_Bulk", "SE_Parent", "SE_Interaction",
                      "Z_intercept", "Z_Bulk", "Z_Parent", "Z_Interaction")) -> GLMdata

GLMdata %>% pivot_wider(names_from = label, values_from = summary) %>% 
  pivot_longer(c(-CHROM, -POS), names_to = c("Param", "Effect"), names_sep = "_") -> GLM_PivotW

saveRDS(GLM_PivotW, file = "Zeocin_ouput.rds")
