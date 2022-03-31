# function (file, highBulk = "HIGH", lowBulk = "LOW", 
#           chromList = NULL, sep = ",") 
# {
  
  # SNPset <- readr::read_delim(file = file, delim = sep, col_names = TRUE, 
  #                             col_types = readr::cols(.default = readr::col_guess(), 
  #                                                     CHROM = "c", POS = "i"))
SNPset <- SelectedDF

  if (!"CHROM" %in% names(SNPset)) {
    stop("No 'CHROM' coloumn found.")
  }
  if (!"POS" %in% names(SNPset)) {
    stop("No 'POS' coloumn found.")
  }

##############################
# ERROR HERE #################
##############################

highBulk <- "NZ.SelectedA.sort"
lowBulk <- "NZ.SelectedC.sort"

  if (!paste0("AD_REF.", highBulk) %in% names(SNPset)) {
    stop("No High Bulk AD_REF coloumn found. Column should be named 'AD_REF.highBulkName'.")
  }

  if (!paste0("AD_REF.", lowBulk) %in% names(SNPset)) {
    stop("No Low Bulk AD_REF coloumn found. Column should be named 'AD_REF.lowBulkName'.")
  }

  if (!paste0("AD_ALT.", highBulk) %in% names(SNPset)) {
    stop("No High Bulk AD_REF coloumn found. Column should be named 'AD_ALT.highBulkName'.")
  }
  if (!paste0("AD_ALT.", lowBulk) %in% names(SNPset)) {
    stop("No Low Bulk AD_ALT coloumn found. Column should be named 'AD_ALT.lowBulkName'.")
  }

##############################
  if (!is.null(chromList)) {
    message("Removing the following chromosomes: ", 
            paste(unique(SNPset$CHROM)[!unique(SNPset$CHROM) %in% 
                                         chromList], collapse = ", "))
    SNPset <- SNPset[SNPset$CHROM %in% chromList, ]
  }
##############################

  SNPset$CHROM <- factor(SNPset$CHROM, levels = gtools::mixedsort(unique(SNPset$CHROM)))
  message("Renaming the following columns: ", paste(colnames(SNPset)[!colnames(SNPset) %in% 
                                                                       c("CHROM", "POS", "REF", "ALT")][grep(highBulk, 
                                                                       x = colnames(SNPset)[!colnames(SNPset) %in% c("CHROM","POS", "REF", "ALT")])], collapse = ", "))
  colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", 
                                            "POS", "REF", "ALT")] <- gsub(pattern = highBulk, 
                                                                          replacement = "HIGH", x = colnames(SNPset)[!colnames(SNPset) %in% 
                                                                                                                       c("CHROM", "POS", "REF", "ALT")])
  message("Renaming the following columns: ", paste(colnames(SNPset)[!colnames(SNPset) %in% 
                                                                       c("CHROM", "POS", "REF", "ALT")][grep(lowBulk, 
                                                                       x = colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")])], collapse = ", "))
  colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", 
                                            "POS", "REF", "ALT")] <- gsub(pattern = lowBulk, 
                                                                          replacement = "LOW", x = colnames(SNPset)[!colnames(SNPset) %in% 
                                                                                                                      c("CHROM", "POS", "REF", "ALT")])
###############################
    SNPset.final <- SNPset %>% dplyr::mutate(DP.HIGH = AD_REF.HIGH + 
                                       AD_ALT.HIGH, DP.LOW = AD_REF.LOW + AD_ALT.LOW, SNPindex.HIGH = AD_ALT.HIGH/DP.HIGH, 
                                     SNPindex.LOW = AD_ALT.LOW/DP.LOW, REF_FRQ = (AD_REF.HIGH + 
                                                                                    AD_REF.LOW)/(DP.HIGH + DP.LOW), deltaSNP = SNPindex.HIGH - 
                                       SNPindex.LOW) %>% dplyr::select(-dplyr::contains("HIGH"), 
                                                                       -dplyr::contains("LOW"), -dplyr::one_of("deltaSNP", 
                                                                                                               "REF_FRQ"), dplyr::matches("AD.*.LOW"), 
                                                                       dplyr::contains("LOW"), dplyr::matches("AD.*.HIGH"), 
                                                                       dplyr::contains("HIGH"), dplyr::everything())
#   return(as.data.frame(SNPset))
# }