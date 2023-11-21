#####
# PROGRAM CODE FOR ARTICLE: SCOTHEART Radiomics data
# (c) Márton Kolossváry, 2023

# DESCRIPTION: Create patient-based data
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
which_type <- "ALL" # "ALL" / "NCP" / "CP"
which_d    <- "EP"
funct      <- mean # Function used for aggregating to patient level data (mean, max)

load(paste0(folder_cache, "Plaque_", which_type, "_", which_d, ".RData"))
p <- readxl::read_xlsx(paste0(folder_data, "SCOT-HEART Per-Patient Masterlist.xlsx"))
setDT(p)
p$subjectno <- as.character(p$subjectno)
p_add <- readxl::read_xlsx(paste0(folder_data, "8y_outcomes_merged_expertreaderplaque.xlsx"))
setDT(p_add)
p_add$Case <- as.character(p_add$Case)
p <- merge(p, p_add[, (c("Case", "MI_cardiac_death_observed_last", "time_MI_cardiac_death_observed_last")), with = FALSE],
           by.x = "subjectno", by.y = "Case", all.x =TRUE, all.y = FALSE)
p$time_MI_cardiac_death_observed_last <- p$time_MI_cardiac_death_observed_last/365

## Create final database -----
l <- merge(p, d_final, by.x = "subjectno", by.y = "patient_ID", all.x = FALSE, all.y = TRUE)
names_eigen <- colnames(l)[(which(colnames(l) == "plaque_ID")+1) : (which(colnames(l) == "Group_name")-1)]
rm(list = c("d_final"))


# SUMMARIZE EIGEN FEATURE ON A PATIENT LEVEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l[, (paste0(names_eigen, "_norm")) := lapply(.SD, function(x) x * volume_1_orig ), .SDcols = names_eigen]
p_eigen <- l[, lapply(.SD, funct), by = "subjectno", .SDcols = paste0(names_eigen, "_norm")]
p_eigen <- cbind(p_eigen[,1], as.data.table(scale(p_eigen[, 2:dim(p_eigen)[2]])))

p <- merge(p, p_eigen, by = "subjectno", all.x = FALSE, all.y = TRUE)
rm(list = c("p_eigen"))


## Format covariates =====
p$risk <- p$assignscore
p$cac  <- log2(p$cacscore +1)
p$obs  <- p$obstructive

p$pb      <- log2(p$Total.TP.burden+1)
p$ncppb   <- log2(p$Total.NCP.burden+1)
p$cppb    <- log2(p$Total.CP.burden+1)
p$lancppb <- log2(p$Total.LD.NCP.burden+1)

p$hasCAD <- p$pb>0


## Save aggregated data =====
save(list = c("p", "l"), file = paste0(folder_cache, "FINAL_", which_type, "_", which_d, ".RData"), envir = .GlobalEnv)
rm(list = ls()); gc(); dev.off()
## ============================================================================
