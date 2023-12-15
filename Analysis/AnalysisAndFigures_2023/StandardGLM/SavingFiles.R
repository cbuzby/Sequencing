#Save all of the data that I've run on AllSequencingRuns.rmd

setwd("C:/Users/cassa/OneDrive/Documents/GitHub/Sequencing/Analysis/AnalysisAndFigures_2023/StandardGLM/")
getwd()
HVYTYDRX2_Fluc_CSS8_glm_prep
HVYTYDRX2_Fluc_CSS1_glm_prep
HVYTYDRX2_CuSO4_CSS8_glm_prep

HJ5HKDRX3_CSS8_glm_prep

HNGLVDRXY_CSS1_glm_prep
xHNGLVDRXY_CSS1_glm_prep

saveRDS(HVYTYDRX2_Fluc_CSS8_glm_prep, file = "HVYTYDRX2_Fluc_CSS8_glm_prep.rds")
saveRDS(HVYTYDRX2_Fluc_CSS1_glm_prep, file = "HVYTYDRX2_Fluc_CSS1_glm_prep.rds")
saveRDS(HVYTYDRX2_CuSO4_CSS8_glm_prep, file = "HVYTYDRX2_CuSO4_CSS8_glm_prep.rds")

saveRDS(HJ5HKDRX3_CSS8_glm_prep, file = "HJ5HKDRX3_CSS8_glm_prep.rds")
saveRDS(HNGLVDRXY_CSS1_glm_prep, file = "HNGLVDRXY_CSS1_glm_prep.rds")
saveRDS(xHNGLVDRXY_CSS1_glm_prep, file = "xHNGLVDRXY_CSS1_glm_prep.rds")

