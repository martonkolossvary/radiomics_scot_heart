#####
# PROGRAM CODE FOR ARTICLE: SCOTHEART Radiomics data
# (c) Márton Kolossváry, 2023

# DESCRIPTION: Normalize and filter radiomic parameters to ones with >0 variance for distance matrix calcualtions
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
which_d <- "EP"
load(paste0(folder_cache, "DATA_", which_type, "_", which_d, ".RData"))

## Model characteristics ------
first_rad  <- 17 #First radiomic parameter column number
# =============================================================================


# NORMALIZE AND REMOVE 0 VARIANCE PARAMETERS %%%%%
## Normalize =====
d[, (dict_rad$term) := lapply(.SD, function(x) {tryCatch(DescTools::Winsorize(x, na.rm = TRUE),
                                                    error = function(e) x)}), .SDcols = dict_rad$term]
d[, (dict_rad$term) := lapply(.SD, function(x) {tryCatch(scale(x),
                                                    error = function(e) x)}), .SDcols = dict_rad$term]

## Clear Inf values =====
for (j in 1:ncol(d)) set(d, which(is.infinite(d[[j]]) | is.nan(d[[j]]) | is.null(d)), j, NA)

## Remove 0 variance columns =====
identical_cols <- sapply(d, function(x) length(unique(x)) == 1)
na_cols        <- sapply(d, function(x) all(is.na(x)))

remove <- union(names(identical_cols)[identical_cols], names(na_cols)[na_cols])

d[, (remove) := NULL]
dict_rad <- dict_rad[!term %in% remove]
# =============================================================================

# SAVE DATA %%%%%
## Save to cache =====
save(list = c("d", "dict_rad"),
     file = paste0(folder_cache, "DATA_", which_type, "_", which_d, "_distmat.RData"), envir = .GlobalEnv)
rm(list = ls()); gc()
# =============================================================================
