#####
# PROGRAM CODE FOR ARTICLE: SCOTHEART Radiomics data
# (c) Márton Kolossváry, 2023

# DESCRIPTION: Correlation plots
#####

# INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load packages =====
library(data.table); library(ggplot2)

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
which_outcome <- 1
limit_CAD <- TRUE

load(paste0(folder_cache, "FINAL_ALL_", which_d, ".RData"))
p_eigen <- colnames(p)[grep("_norm", colnames(p),fixed = TRUE)]
p_eigen <- p_eigen[-length(p_eigen)]
eigen_all <- p[, ..p_eigen]
colnames(eigen_all) <- paste0(colnames(eigen_all), "_ALL")

load(paste0(folder_cache, "FINAL_NCP_", which_d, ".RData"))
p_eigen <- colnames(p)[grep("_norm", colnames(p),fixed = TRUE)]
p_eigen <- p_eigen[-length(p_eigen)]
eigen_ncp <- p[, ..p_eigen]
colnames(eigen_ncp) <- paste0(colnames(eigen_ncp), "_NCP")

load(paste0(folder_cache, "FINAL_CP_", which_d, ".RData"))
p_eigen <- colnames(p)[grep("_norm", colnames(p),fixed = TRUE)]
p_eigen <- p_eigen[-length(p_eigen)]
eigen_cp <- p[, ..p_eigen]
colnames(eigen_cp) <- paste0(colnames(eigen_cp), "_CP")

eigen <- cbind(eigen_all, eigen_ncp, eigen_cp) # ALL Radiomic parameters
p <- cbind(p, eigen)

if(limit_CAD) {
 p <- p[pb >0, ]
}

vars <- c(c("risk", "cac", "obs", "Total.NCP.vol..mm3.", "Total.CP.vol..mm3.", "Total.LD.NCP.volume..mm3."), colnames(eigen))
vars <- c(c("ager", "gendern", "bmir", "hyp", "dm", "hyperlipid", "smoking", "risk", "cac", "obs", "Total.NCP.vol..mm3.", "Total.CP.vol..mm3.", "Total.LD.NCP.volume..mm3."), colnames(eigen))
outcomes <- c("MI_cardiac_death_observed_last")
times    <- c("time_MI_cardiac_death_observed_last")
# =============================================================================


# PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggcorrplot)
### Functions
nice_p <- function(p) {
 ifelse(p<0.001, "p<0.001", ifelse(p<0.01, "p<0.01", sprintf("p=%.2f", round(p, 2))))
}
pretty_table <- function(d, row_names = v[1:dim(d)[1]]){
 d$Estimate <- sprintf("%.2f", round((10^d$estimate -1)*100, 2))
 d$Confidence_interval <- paste0("[", sprintf("%.2f", round((10^d$conf.low -1)*100, 2)), "; ", sprintf("%.2f", round((10^d$conf.high -1)*100, 2)), "]")
 d$p <- format.pval(d$p.value, digits = 1, eps = 0.001, na.form = "NA", nsmall = 3)
 d$Variable <- row_names
 d
}


## Plot correlation -----
p_data <- p[, ..vars]
vars <- c(c("ager", "gendern", "bmir", "hyp", "dm", "hyperlipid", "smoking", "risk", "cac", "obs",
            "Total.NCP.vol..mm3.", "Total.CP.vol..mm3.", "Total.LD.NCP.volume..mm3.",
            "ncppb", "cppb", "lancppb",
            "brown_norm_ALL", "black_norm_ALL", "greenyellow_norm_NCP", "pink_norm_NCP",
            "magenta_norm_NCP", "brown_norm_NCP", "yellow_norm_NCP", "brown_norm_CP"))

vars <- c(c("risk", "cac", "obs",
            "Total.NCP.vol..mm3.", "Total.CP.vol..mm3.", "Total.LD.NCP.volume..mm3.",
            "ncppb", "cppb", "lancppb",
            "brown_norm_ALL", "black_norm_ALL", "greenyellow_norm_NCP", "pink_norm_NCP",
            "magenta_norm_NCP", "brown_norm_NCP", "yellow_norm_NCP", "brown_norm_CP"))

d_cor <- p[, ..vars]
cor_mat <- cor(d_cor, use = "pairwise.complete.obs", method = "spearman") # pearson / spearman
colnames(cor_mat) = row.names(cor_mat) <- c("Age", "Sex", "BMI", "Hypertension","Diabetes", "Hyperlipidaemia", "Smoking", "ASSIGN risk score", "Log2 Agatston-score", "Obstructive CAD",
                                            "NCP  volume", "CP  volume", "LANCP volume", "log2 NC-PB", "log2 C-PB", "log2 LANC-PB",
                                            "Eigen-All 3", "Eigen-All 5", "Eigen-NCP 3", "Eigen-NCP 4", "Eigen-NCP 5", "Eigen-NCP 6", "Eigen-NCP 7", "Eigen-CP 1")

colnames(cor_mat) = row.names(cor_mat) <- c("ASSIGN risk score", "Log2 Agatston-score", "Obstructive disease", "NCP  volume", "CP  volume", "LANCP volume", "log2 NC-PB", "log2 C-PB", "log2 LANC-PB",
                                            "Eigen-All 3", "Eigen-All 5", "Eigen-NCP 3", "Eigen-NCP 4", "Eigen-NCP 5", "Eigen-NCP 6", "Eigen-NCP 7", "Eigen-CP 1")
cor_mat[is.na(cor_mat)] <- 0
## Plot -----
cor_p <- ggcorrplot::ggcorrplot(cor_mat, hc.order = TRUE, lab = TRUE,
                                title = "Correlation between clinical factors, volumetric and significant eigen radiomic features", # "Correlation between proteins, biomarkers and NCP plaque changes"
                                legend.title =  "Spearman correlation", # "Pearson correlation"
                                colors = c("#6D9EC1", "white", "#E46726"), lab_size = 3)
ggplot2::ggsave(path = folder_img, filename = "Correlation_clin_plaque.pdf", plot = cor_p, device = "pdf", width = 15, height = 12.5, units = "in") # 12/10 18/15 15/12.5