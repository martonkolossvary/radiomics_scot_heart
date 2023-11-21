#####
# PROGRAM CODE FOR ARTICLE: SCOTHEART Radiomics data
# (c) Márton Kolossváry, 2023

# DESCRIPTION: Create lesion level radiomics
#####

# INITIALIZE ====================
## Load packages -----
library(RIA); library(foreach); library(doParallel)

## Get files -----
folder_nii  <- "" 
folder_csv  <- "" 
files_long  <- list.files(folder_nii, pattern = ".nii", full.names = TRUE)
files_short <- list.files(folder_nii, pattern = ".nii", full.names = FALSE)
which_type <- "ALL" # "ALL" / "NCP" / "CP"
ncores <- 20

## Remove already done -----
keep <- NULL
for(i in 1:length(files_short)) {
 if(!file.exists(paste0(folder_csv, "L_", which_type, "_", substr(files_short[i], 1, nchar(files_short[i])-4), ".csv"))) {
  keep <- c(keep, i)
 }
}
files_long  <- files_long[keep]
files_short <- files_short[keep]

# CALCULATE PARALLEL ====================
## Initiate threads -----
doParallel::registerDoParallel(ncores)

data_out_paral <- foreach(i = sample.int(length(files_short), n = length(files_short), replace = FALSE), .combine="rbind", .inorder=TRUE,
                          .packages=c('RIA'), .errorhandling = c("remove"), .verbose=TRUE) %dopar% {
                           
                           if(which_type == "ALL") {
                            if(!file.exists(paste0(folder_csv, "L_ALL_", substr(files_short[i], 1, nchar(files_short[i])-4), ".csv"))) {
                             d <- RIA::load_nifti(files_long[i], verbose_in = FALSE)
                             if(d$log$orig_vol_mm != 0){
                              d <- RIA::radiomics_all(d, verbose_in = FALSE)
                              RIA::save_RIA(d, save_to = folder_csv, save_name = paste0("L_ALL_", substr(files_short[i], 1, nchar(files_short[i])-4)), group_name = files_short[i])
                             }
                             rm(d); gc()
                            }
                           }
                           
                           if(which_type == "NCP") {
                            if(!file.exists(paste0(folder_csv, "L_NCP_", substr(files_short[i], 1, nchar(files_short[i])-4), ".csv"))) {
                             d <- RIA::load_nifti(files_long[i], mask_filename = files_long[i], keep_mask_values = "<=350", verbose_in = FALSE)
                             if(d$log$orig_vol_mm != 0){
                              d <- RIA::radiomics_all(d, verbose_in = FALSE)
                              RIA::save_RIA(d, save_to = folder_csv, save_name = paste0("L_NCP_", substr(files_short[i], 1, nchar(files_short[i])-4)), group_name = files_short[i])
                             }
                             rm(d); gc()
                            }
                           }
                           
                           if(which_type == "CP") { 
                            if(!file.exists(paste0(folder_csv, "L_CP_", substr(files_short[i], 1, nchar(files_short[i])-4), ".csv"))) {
                             d <- RIA::load_nifti(files_long[i], mask_filename = files_long[i], keep_mask_values = ">350", verbose_in = FALSE)
                             if(d$log$orig_vol_mm != 0){
                              d <- RIA::radiomics_all(d, verbose_in = FALSE)
                              RIA::save_RIA(d, save_to = folder_csv, save_name = paste0("L_CP_", substr(files_short[i], 1, nchar(files_short[i])-4)), group_name = files_short[i])
                             }
                             rm(d); gc()
                            }
                           }
                          }
