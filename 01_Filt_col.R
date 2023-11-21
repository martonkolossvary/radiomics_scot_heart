#####
# PROGRAM CODE FOR ARTICLE: SCOTHEART Radiomics data
# (c) Márton Kolossváry, 2023

# DESCRIPTION: Create databases with only previously published radiomic features
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
d <- data.table::fread(paste0(folder_data, "FINAL/LESION_", which_type, ".csv"))

## Model characteristics ------
first_rad  <- 17 #First radiomic parameter column number
type       <- 2 #Type of discretizations
bins       <- c(8, 16, 32) #Number of bins
dist       <- 1 #Distances for GLCM
geom_keep  <- rep(TRUE, 12) #Needed geometric statistics
# =============================================================================


# CREATE PATIENT IDs ====================
patient_ID <- substr(d$Group_name, 1, regexpr("__", d$Group_name)-1)
plaque_ID  <- substr(d$Group_name, regexpr("__", d$Group_name)+2, nchar(d$Group_name)-4)

d[, V1:=NULL]; d[, V1:=NULL]
d <- cbind(patient_ID, plaque_ID, d)

# RADIOMIC FEATURE SETS =====
rad_params <- colnames(d)[first_rad:dim(d)[2]]
rad_params <- janitor::make_clean_names(rad_params)
colnames(d)[first_rad:dim(d)[2]] <- rad_params

## First order -----
fo <- 1:44

## GLCM ------
glcm_start <- fo[length(fo)]+1
glcm       <- NULL
for (i in 1:(type*length(bins)*dist)) {
 glcm_codes_3 <- seq((glcm_start+240*(i-1)), (glcm_start+240*(i-1)+188), 3)
 glcm_codes_1 <- (glcm_start+240*(i-1)+189):(glcm_start+240*i-1)
 glcm         <- c(glcm, glcm_codes_3, glcm_codes_1)
}

## GLRLM -----
glrlm_start <- glcm[length(glcm)]+1
glrlm       <- glrlm_start:(glrlm_start + 11*type*length(bins)*dist-1)

## Geometric -----
geom_start <- glrlm[length(glrlm)]+1
bins_all   <- rep(bins, type)
geom_boo   <- geom_keep #Calculate on original image also
for (i in 1:length(bins_all)) {
 for(j in 1:length(geom_keep)) {
  geom_boo <- c(geom_boo, rep(geom_keep[j], bins_all[i]) )
 }
}
geom <- geom_start:(geom_start + length(geom_boo) -1)
geom <- geom[geom_boo]

keep_rad <- rep(FALSE, length(rad_params))
keep_rad[c(fo, glcm, glrlm, geom)] <- TRUE
# =============================================================================


# CREATE DATA TABLE ACCORDING TO RADIOMICS FEATURES ======
fo_c    <- rep("First_order", length(fo))
glcm_c  <- rep("GLCM", length(glcm[1]:(glrlm[1]-1)))
glrlm_c <- rep("GLRLM", length(glrlm[1]:(geom[1]-1)))
geom_c  <- rep("Geometrical", length(geom[1]:(geom[length(geom)])))

rad_class <- c(fo_c, glcm_c, glrlm_c, geom_c)

## Bin type elimination -----
ep_rgx <- "_ep_"
es_rgx <- "_es_"

ep <- grepl(ep_rgx, rad_params, fixed = TRUE)
es <- grepl(es_rgx, rad_params, fixed = TRUE)

## Bin number elimination -----
b8_rgx  <- "_b8_";  b8_rgx_geom  <- "_8$";  b8_rgx_glcm  <- "_b8$"
b16_rgx <- "_b16_"; b16_rgx_geom  <- "_16$"; b16_rgx_glcm  <- "_b16$"
b32_rgx <- "_b32_"; b32_rgx_geom  <- "_32$"; b32_rgx_glcm  <- "_b32$"

b8  <- grepl(b8_rgx, rad_params, fixed = TRUE) | grepl(b8_rgx_geom, rad_params, fixed = FALSE) | grepl(b8_rgx_glcm, rad_params, fixed = FALSE)
b16 <- grepl(b16_rgx, rad_params, fixed = TRUE) | grepl(b16_rgx_geom, rad_params, fixed = FALSE) | grepl(b16_rgx_glcm, rad_params, fixed = FALSE)
b32 <- grepl(b32_rgx, rad_params, fixed = TRUE) | grepl(b32_rgx_geom, rad_params, fixed = FALSE) | grepl(b32_rgx_glcm, rad_params, fixed = FALSE)
# =============================================================================


# CREATE AND SAVE OUTPUTS ======
dict_rad <- as.data.table(cbind(rad_params, rad_class, keep_rad, es, ep, b8, b16, b32))
colnames(dict_rad)[1] <- "term"
change_boo <- colnames(dict_rad)[3:8]
dict_rad[, (change_boo) := lapply(.SD, as.logical), .SDcols = change_boo]


## Save Radiomics data -----
d_all <- copy(d); dict_rad_all <- copy(dict_rad)
d_rad <- d[, first_rad:dim(d)[2]]

KEEP <- c("keep_rad", "ep")
KEEP  <- apply(dict_rad_all[, ..KEEP], 1, sum, na.rm = T) == length(KEEP) | apply(dict_rad_all[, grep("[0-9]", colnames(dict_rad_all)), with = FALSE], 1, sum) == 0
dict_rad <- dict_rad_all[KEEP, ]
d <- cbind(d_all[, 1:(first_rad-1)], d_rad[, ..KEEP])
save(list = c("d", "dict_rad"),
     file = paste0(folder_cache, "DATA_", which_type, "_EP.RData"), envir = .GlobalEnv)

rm(list = ls()); gc()
# =============================================================================
