#####
# PROGRAM CODE FOR ARTICLE: SCOTHEART Radiomics data
# (c) Márton Kolossváry, 2022

# DESCRIPTION: Cox analysis
#####

# INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load packages =====
library(data.table)

## File and folder locations =====
set.seed(42)
folder_wd    <- ""
folder_data  <- paste0(folder_wd, "DATA/")
folder_stat  <- paste0(folder_wd, "STAT/")
folder_cache <- paste0(folder_wd, "CACHE/")
folder_img   <- paste0(folder_wd, "IMAGES/")
folder_src   <- paste0(folder_wd, "SCRIPTS/")

## Load data =====
which_d    <- "EP"
prefix <- "COX_Vol_FINAL_"

d_ALL <- fread(paste0(folder_stat, "RAW/", prefix, "ALL_", which_d, ".csv"))
d_NCP <- fread(paste0(folder_stat, "RAW/", prefix, "NCP_", which_d, ".csv"))
d_CP <- fread(paste0(folder_stat, "RAW/", prefix, "CP_", which_d, ".csv"))


# CREATE OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for (outcome in unique(d_ALL$Outcome)) {
 out_ALL <- data.table(EigenRadiomics = unique(d_ALL$EigenRadiomics))
 out_NCP <- data.table(EigenRadiomics = unique(d_NCP$EigenRadiomics))
 out_CP  <- data.table(EigenRadiomics = unique(d_CP$EigenRadiomics))
 
 for (model in  unique(d_ALL$Covariates)) {
  d_i_ALL <- d_ALL[Outcome == outcome & Covariates == model]
  out_ALL <- cbind(out_ALL, data.table(HR = sprintf("%.2f", round(d_i_ALL$HR, 2)),
                                       CI = paste0("[", sprintf("%.2f", round(d_i_ALL$ConfLow, 2)), "; ", sprintf("%.2f", round(d_i_ALL$ConfHigh, 2)), "]"),
                                       p  = sprintf("%.4f", round(d_i_ALL$p, 4))))
  d_i_NCP <- d_NCP[Outcome == outcome & Covariates == model]
  out_NCP <- cbind(out_NCP, data.table(HR = sprintf("%.2f", round(d_i_NCP$HR, 2)),
                                       CI = paste0("[", sprintf("%.2f", round(d_i_NCP$ConfLow, 2)), "; ", sprintf("%.2f", round(d_i_NCP$ConfHigh, 2)), "]"),
                                       p  = sprintf("%.4f", round(d_i_NCP$p, 4))))
  d_i_CP <- d_CP[Outcome == outcome & Covariates == model]
  out_CP <- cbind(out_CP, data.table(HR = sprintf("%.2f", round(d_i_CP$HR, 2)),
                                     CI = paste0("[", sprintf("%.2f", round(d_i_CP$ConfLow, 2)), "; ", sprintf("%.2f", round(d_i_CP$ConfHigh, 2)), "]"),
                                     p  = sprintf("%.4f", round(d_i_CP$p, 4))))
 }
 info <- data.table(EigenRadiomics = c("", paste0("Outcome: ", outcome), unlist(lapply(1:length(unique(d_ALL$Covariates)), function(x) {paste0("Model-", x, ": ", unique(d_ALL$Covariates)[x])}))))
 out_ALL <- rbindlist(list(out_ALL, info), fill = TRUE, )
 out_NCP <- rbindlist(list(out_NCP, info), fill = TRUE)
 out_CP <- rbindlist(list(out_CP, info), fill = TRUE)
 
 write.csv(out_ALL, paste0(folder_stat, prefix, which_d, "_", outcome, "_ALL.csv"), na = "", row.names = FALSE)
 write.csv(out_NCP, paste0(folder_stat, prefix, which_d, "_", outcome, "_NCP.csv"), na = "", row.names = FALSE)
 write.csv(out_CP, paste0(folder_stat, prefix, which_d, "_", outcome, "_CP.csv"), na = "", row.names = FALSE)
}

rm(list = ls()); gc(); dev.off()
