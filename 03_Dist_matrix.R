#####
# PROGRAM CODE FOR ARTICLE: SCOTHEART Radiomics data
# (c) Márton Kolossváry, 2023

# DESCRIPTION: Create distance matrix for WGCNA considering using LMM with patient level random intercept
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
load(paste0(folder_cache, "DATA_", which_type, "_", which_d, "_distmat.RData"))

## Model characteristics ------
first_rad  <- 8 #First radiomic parameter column number
DISTANCE   <- "LMM" # Distance metric between pairs of radiomic parameters: euclidean, pearson, spearman, kendall, LMM
patient_ID <- "patient_ID" 
SQUARE     <- FALSE # Square the corresponding distance (i.e. to get R2 from pearson)
n_cores    <- 14 # Number of cores to use for calculation in case of LMM
# =============================================================================


# CREATE DISTANCE MATRIX =====
if(DISTANCE == "euclidean") {
 m_dist <- rdist::pdist(t(d[, first_rad:(dim(d)[2])]), metric = DISTANCE)
 save(m_dist, file = paste0(folder_cache, "dist_mat_", which_type, "_", which_d, ".RData"))
}

if(DISTANCE %in% c("pearson", "spearman", "kendall")) {
 m_dist <- cor(d[, first_rad:(dim(d)[2])], method = DISTANCE)
 m_dist[is.na(m_dist)] <- 0
 m_dist <- abs(m_dist)
 if(SQUARE) {m_dist <- m_dist^2}
 save(m_dist, file = paste0(folder_cache, "dist_mat_", which_type, "_", which_d, ".RData"))
}

if(DISTANCE == "LMM") {
 #https://www.pipinghotdata.com/posts/2021-10-11-estimating-correlations-adjusted-for-group-membership/
 #https://stats.stackexchange.com/questions/111150/calculating-r2-in-mixed-models-using-nakagawa-schielzeths-2013-r2glmm-me
 
 library(foreach); library(doParallel); library(bigstatsr)
 doParallel::registerDoParallel(n_cores)
 t_start <- Sys.time()
 print(paste("Starting calcualtion of distance matrix at:", Sys.time()))
 
 ## Initiate matrix ----
 d_rad_lmm <- d[, first_rad:(dim(d)[2])]
 FBM_dist       <- FBM(nrow = dim(d_rad_lmm)[2], ncol = dim(d_rad_lmm)[2], type = "double", init = 0)
 diag(FBM_dist) <- 1
 
 d_rad_lmm$GROUP_VAR <- d[[patient_ID]]
 
 ## Start calculations -----
 junk <- foreach(i = 1:(dim(d_rad_lmm)[2]-2), .combine='c') %:%
  foreach(j = (i+1):(dim(d_rad_lmm)[2]-1), .combine='c') %dopar% {
   model        <- lmerTest::lmer(paste0(colnames(d_rad_lmm)[i], " ~ ", colnames(d_rad_lmm)[j], " + (1|GROUP_VAR)"), data = d_rad_lmm)
   
   if(performance::r2_nakagawa(model, by_group = TRUE)$R2[1] < 0) {
    FBM_dist[i, j] <- 0
   } else {
    if(SQUARE) {
     FBM_dist[i, j] <- performance::r2_nakagawa(model, by_group = TRUE)$R2[1]
    } else {
     FBM_dist[i, j] <- sqrt(performance::r2_nakagawa(model, by_group = TRUE)$R2[1])
    }
   }
   NULL
  }; rm(junk)
 
 m_dist <- FBM_dist[1:dim(FBM_dist)[1], 1:dim(FBM_dist)[2]]
 rm(FBM_dist)
 
 m_dist[lower.tri(m_dist, diag = FALSE)] <- 0
 t_m_dist <- t(m_dist); diag(t_m_dist) <- 0
 m_dist <- m_dist + t_m_dist; rm(t_m_dist)
 m_dist[m_dist<0] <- 0
 m_dist[is.na(m_dist) | is.nan(m_dist) | is.infinite(m_dist)] <- 0
 
 save(m_dist, file = paste0(folder_cache, "dist_mat_", which_type, "_", which_d, ".RData"))
 print(paste("Time taken to build Distance matrix:", Sys.time() - t_start, "hours"))
}
# =============================================================================
rm(list = ls()); gc(); dev.off()