#####
# PROGRAM CODE FOR ARTICLE: SCOTHEART Radiomics data
# (c) Márton Kolossváry, 2023

# DESCRIPTION: Multivariate Cox models
#####

# INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load packages =====
library(data.table); library(survival); library(survminer); library(randomForestSRC)

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
prefix <- "COX_logPB_FINAL"
which_outcome <- 1

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
# write.csv(l[order(l$magenta_norm)[c(1:5, 4038:4042)], c("subjectno", "id", "plaque_ID", "Total.plaque.vol..mm3.", "magenta_norm")], "/Users/mjk2/Downloads/Magenta_NCP.csv", row.names = F)


load(paste0(folder_cache, "FINAL_CP_", which_d, ".RData"))
p_eigen <- colnames(p)[grep("_norm", colnames(p),fixed = TRUE)]
p_eigen <- p_eigen[-length(p_eigen)]
eigen_cp <- p[, ..p_eigen]
colnames(eigen_cp) <- paste0(colnames(eigen_cp), "_CP")
# write.csv(l[order(l$brown_norm)[c(1:5, 4038:4042)], c("subjectno", "id", "plaque_ID", "Total.plaque.vol..mm3.", "brown_norm")], "/Users/mjk2/Downloads/Brown_CP.csv", row.names = F)

eigen <- cbind(eigen_all, eigen_ncp, eigen_cp) # ALL Radiomic parameters
only_sig <- list(c("brown_norm_ALL", "black_norm_ALL", "greenyellow_norm_NCP", "pink_norm_NCP",
                   "magenta_norm_NCP", "brown_norm_NCP", "yellow_norm_NCP", "brown_norm_CP"))[[which_outcome]] # Plaque burden 8 year

eigen <- eigen[, ..only_sig, with = FALSE]
p <- cbind(p, eigen)

p$Total.NCP.vol..mm3. <- p$Total.NCP.vol..mm3./100
p$Total.CP.vol..mm3. <- p$Total.CP.vol..mm3./100
p$Total.LD.NCP.volume..mm3. <- p$Total.LD.NCP.volume..mm3./100

# covars <- c("Total.NCP.vol..mm3.", "Total.CP.vol..mm3.", "Total.LD.NCP.volume..mm3.", "risk", "cac", "obs")
covars <- c("ncppb", "cppb", "lancppb", "risk", "cac", "obs")
vars <- c(covars, colnames(eigen))

outcomes <- c("MI_cardiac_death_observed_last")
times    <- c("time_MI_cardiac_death_observed_last")
# =============================================================================

# RUN MUTLIVARIATE COX ====================
mod <- survival::coxph(
 formula(paste0("survival::Surv(time =", times[which_outcome], ", event=", outcomes[which_outcome], ") ~ ", paste0(vars, collapse = "+"))), data = p)
mod0 <- survival::coxph(
 formula(paste0("survival::Surv(time =", times[which_outcome], ", event=", outcomes[which_outcome], ") ~ ", paste0(covars, collapse = "+"))), data = p)
mod00 <- survival::coxph(
 formula(paste0("survival::Surv(time =", times[which_outcome], ", event=", outcomes[which_outcome], ") ~ ", paste0(covars[4:6], collapse = "+"))), data = p)

lmtest::lrtest(mod, mod0)
lmtest::lrtest(mod0, mod00)

which_mod <- mod
summ <- broomExtra::tidy(which_mod, exponentiate = TRUE, conf.int = TRUE)[, c(1, 2, 6, 7, 5)]
summ$estimate <- sprintf("%.2f", round(summ$estimate, 2))
summ$CI   <- paste0("[", sprintf("%.2f", round(summ$conf.low, 2)), "; ", sprintf("%.2f", round(summ$conf.high, 2)), "]")
summ$p.value <- sprintf("%.4f", round(summ$p.value, 4))
summ$conf.low <- NULL; summ$conf.high <- NULL;
summ <- summ[, c(1,2,4,3)]
#write.csv(summ, paste0(folder_stat, prefix, which_d, "_", outcomes[which_outcome], "_multivariate_FDR.csv"), na = "", row.names = FALSE)