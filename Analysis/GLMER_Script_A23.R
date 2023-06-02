#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Load packages
library(tidyr)
library(reshape2)
library(dplyr)
library(lme4)

#Load in function
glm_mixed_model_Fluc6_div <- function(D_OI_B_Oak, D_OI_B_Wine, D_WI_B_Oak,D_WI_B_Wine,
                                      D_OI_C_Oak,D_OI_C_Wine,D_OI_D_Oak, D_OI_D_Wine,
                                      F_OI_B_Oak, F_OI_B_Wine, F_OI_C_Oak,F_OI_C_Wine,
                                      F_OI_D_Oak,F_OI_D_Wine,F_WI_C_Oak, F_WI_C_Wine,
                                      F_WI_D_Oak, F_WI_D_Wine){
  
  combineddata <- data.frame(Bulk = factor(c("L", "L","L", "L", "L", "L","L", "L", 
                                             "H", "H","H", "H", "H", "H","H", "H", 
                                             "H", "H")),
                             Parent = factor(c("O", "O", "W", "W", "O", "O", "O", "O",  
                                               "O", "O", "O", "O", "O", "O", "W", "W", 
                                               "W", "W")),
                             Allele = factor(c("O", "W", "O", "W","O", "W", "O", "W", 
                                               "O", "W", "O", "W","O", "W", "O", "W",
                                               "O", "W")),
                             Rep = factor(1:18),
                             Reads = c(D_OI_B_Oak, D_OI_B_Wine, D_WI_B_Oak,D_WI_B_Wine,
                                       D_OI_C_Oak,D_OI_C_Wine,D_OI_D_Oak, D_OI_D_Wine,
                                       F_OI_B_Oak, F_OI_B_Wine, F_OI_C_Oak,F_OI_C_Wine,
                                       F_OI_D_Oak,F_OI_D_Wine,F_WI_C_Oak, F_WI_C_Wine,
                                       F_WI_D_Oak, F_WI_D_Wine))
  
  combineddata %>% group_by(Bulk, Parent) %>% summarize(key = factor(paste(Bulk, Parent, sep = "_"))) %>% merge(combineddata) -> combineddata
  
  b <- glmer(Allele ~ Bulk*Parent+(1+key | Rep), 
             weights = Reads, family = binomial, 
             data = combineddata)
  
  
  return(c(summary(b)$coefficients[1:12]))
  
} 

#Load data
cybr2Data <- readRDS("Fluconazole_1_cybr2.rds")

#Roll data
cybr2Data %>% pivot_longer(c(-CHROM, -POS), names_to = "label") %>% 
  group_by(CHROM, label) %>% arrange(POS) %>% 
  summarize(POS = POS, 
            CHROM = CHROM, 
            SmoothCount = ceiling(frollapply(value, n = 200, FUN = median))) %>% #CHANGE THIS LINE
  na.omit() %>% 
  pivot_wider(names_from = label,values_from = SmoothCount) -> rollData


#Run function
rollData %>% group_by(CHROM, POS) %>% #head(1) %>%
  summarise(summary = glm_mixed_model_Fluc6_div(D_OI_B_Oak = Dilute_OakI_B_Oak, 
                                                D_OI_B_Wine = Dilute_OakI_B_Wine, 
                                                D_WI_B_Oak = Dilute_WineI_B_Oak,
                                                D_WI_B_Wine = Dilute_WineI_B_Wine,
                                                
                                                D_OI_C_Oak = Dilute_WineI_C_Oak,
                                                D_OI_C_Wine = Dilute_WineI_C_Wine,
                                                D_OI_D_Oak = Dilute_WineI_D_Oak, 
                                                D_OI_D_Wine = Dilute_WineI_D_Wine,
                                                
                                                F_OI_B_Oak = Fluconazole_OakI_B_Oak, 
                                                F_OI_B_Wine = Fluconazole_OakI_B_Wine, 
                                                F_OI_C_Oak = Fluconazole_OakI_C_Oak,
                                                F_OI_C_Wine = Fluconazole_OakI_C_Wine,
                                                
                                                F_OI_D_Oak = Fluconazole_OakI_D_Oak,
                                                F_OI_D_Wine = Fluconazole_OakI_D_Wine,
                                                F_WI_C_Oak = Fluconazole_WineI_C_Oak, 
                                                F_WI_C_Wine = Fluconazole_WineI_C_Wine,
                                                
                                                F_WI_D_Oak = Fluconazole_WineI_D_Oak,
                                                F_WI_D_Wine = Fluconazole_WineI_D_Wine),
            label = c("E_intercept", "E_Bulk", "E_Parent", "E_Interaction",
                      "SE_intercept", "SE_Bulk", "SE_Parent", "SE_Interaction",
                      "Z_intercept", "Z_Bulk", "Z_Parent", "Z_Interaction")) -> GLMdata

GLMdata %>% pivot_wider(names_from = label, values_from = summary) %>% 
  pivot_longer(c(-CHROM, -POS), names_to = c("Param", "Effect"), names_sep = "_") -> GLM_PivotW

saveRDS(GLM_PivotW, file = "Fluconazole_ouput.rds")
