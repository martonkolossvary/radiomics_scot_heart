#####
# PROGRAM CODE FOR ARTICLE: SCOTHEART Radiomics data
# (c) Márton Kolossváry, 2022

# DESCRIPTION: ROC and KM curves
#####

# INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load packages =====
library(data.table); library(survival); library(survminer); library(pROC); library(ggplot2)

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

load(paste0(folder_cache, "FINAL_CP_", which_d, ".RData"))
p_eigen <- colnames(p)[grep("_norm", colnames(p),fixed = TRUE)]
p_eigen <- p_eigen[-length(p_eigen)]
eigen_cp <- p[, ..p_eigen]
colnames(eigen_cp) <- paste0(colnames(eigen_cp), "_CP")

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

# p = p[gender == "Male"]
# RUN MUTLIVARIATE COX ====================
mod <- survival::coxph(
 formula(paste0("survival::Surv(time =", times[which_outcome], ", event=", outcomes[which_outcome], ") ~ ", paste0(vars, collapse = "+"))), data = p)
mod0 <- survival::coxph(
 formula(paste0("survival::Surv(time =", times[which_outcome], ", event=", outcomes[which_outcome], ") ~ ", paste0(covars, collapse = "+"))), data = p)
mod00 <- survival::coxph(
 formula(paste0("survival::Surv(time =", times[which_outcome], ", event=", outcomes[which_outcome], ") ~ ", paste0(covars[4:6], collapse = "+"))), data = p)

roc00 <- roc(p[[outcomes]], predict(mod00))
roc0 <- roc(p[[outcomes]], predict(mod0))
roc <- roc(p[[outcomes]], predict(mod))

roc_comp0 <- roc.test(roc00, roc0)
roc_comp <- roc.test(roc0, roc)

# PLOT ROC CURVES =====
p_roc <- ggroc(list(roc00, roc0, roc)) +
 xlab("Specificity") + ylab("Sensitivity") + 
 geom_abline(slope=1, intercept = 1, linetype = "dashed", color = "grey") +
 theme_bw() + theme(legend.position = "bottom") + coord_fixed(ratio = 1.0) +
 labs(title = "ROC of multivariate Cox regression models", color = "Models") +
 scale_color_manual(labels = c("Clinical", "Clinical + Plaque burden", "Clinical + Plaque burden + Eigen radiomics"),
                    values = c("#1f77b4ff", "#ff7f0eff", "#2ca02cff")) +
 annotate("text", x=0.45, y=0.15, label= "AUC: 0.72, 95%CI: 0.66-0.77", color = "#1f77b4ff", size = 6) +
 annotate("text", x=0.45, y=0.10, label= "AUC: 0.73, 95%CI: 0.68-0.78", color = "#ff7f0eff", size = 6) +
 annotate("text", x=0.45, y=0.05, label= "AUC: 0.76, 95%CI: 0.71-0.80", color = "#2ca02cff", size = 6)
 
ggplot2::ggsave(path = folder_img, filename = "ROC_pb.pdf", plot = p_roc, device = "pdf", width = 8, height = 8, units = "in") # 12/10 18/15 15/12.5


p_roc <- ggroc(list(roc00, roc0, roc)) +
 xlab("Specificity") + ylab("Sensitivity") + 
 geom_abline(slope=1, intercept = 1, linetype = "dashed", color = "grey") +
 theme_bw() + theme(legend.position = "bottom") + coord_fixed(ratio = 1.0) +
 labs(title = "ROC of multivariate Cox regression models", color = "Models") +
 scale_color_manual(labels = c("Clinical", "Clinical + Plaque volume", "Clinical + Plaque volume + Eigen radiomics"),
                    values = c("#1f77b4ff", "#ff7f0eff", "#2ca02cff")) +
 annotate("text", x=0.45, y=0.15, label= "AUC: 0.72, 95%CI: 0.66-0.77", color = "#1f77b4ff", size = 6) +
 annotate("text", x=0.45, y=0.10, label= "AUC: 0.71, 95%CI: 0.66-0.77", color = "#ff7f0eff", size = 6) +
 annotate("text", x=0.45, y=0.05, label= "AUC: 0.75, 95%CI: 0.70-0.80", color = "#2ca02cff", size = 6)

ggplot2::ggsave(path = folder_img, filename = "ROC_vol.pdf", plot = p_roc, device = "pdf", width = 8, height = 8, units = "in") # 12/10 18/15 15/12.5


### KM-curves regression -----
library(survival); library(survminer)

roc <- roc(p[[outcomes]], p$ncppb)
cord_ncppb <- coords(roc, "best")
p$ncppb_med <- factor(as.numeric(p$ncppb>cord_ncppb$threshold))

roc <- roc(p[[outcomes]], p$cppb)
cord_cppb <- coords(roc, "best")
p$cppb_med <- factor(as.numeric(p$cppb>cord_cppb$threshold))

roc <- roc(p[[outcomes]], p$lancppb)
cord_lancppb <- coords(roc, "best")
p$lancppb_med <- factor(as.numeric(p$lancppb>cord_lancppb$threshold))

roc <- roc(p[[outcomes]], p$risk)
cord_risk <- coords(roc, "best")
p$risk_med <- factor(as.numeric(p$risk>cord_risk$threshold))

roc <- roc(p[[outcomes]], p$cac)
cord_cac <- coords(roc, "best")
p$cac_med <- factor(as.numeric(p$cac>cord_cac$threshold))

p$obs_med <- factor(as.numeric(abs(p$obs-1)))

roc <- roc(p[[outcomes]], p$magenta_norm_NCP)
cord_magenta_norm_NCP <- coords(roc, "best")
p$magenta_norm_NCP_med <- factor(as.numeric(p$magenta_norm_NCP>cord_magenta_norm_NCP$threshold))

roc <- roc(p[[outcomes]], p$brown_norm_CP)
cord_brown_norm_CP <- coords(roc, "best")
p$brown_norm_CP_med <- factor(as.numeric(p$brown_norm_CP>cord_brown_norm_CP$threshold))
 

fit <- survfit(Surv(p[[times[which_outcome]]], p[[outcomes[which_outcome]]]) ~ risk_med, data = p)
p_1 <- ggsurvplot(fit, pval = TRUE, pval.method = TRUE, risk.table = FALSE, test.for.trend = FALSE, conf.int = FALSE,
                  xlab = "Years", palette = c("#00a1aa", "#aa004c"), ylim = c(0.87,1), pval.coord = c(2, 0.87), pval.method.coord = c(0,0.87),
                  legend.labs = c(paste0(">", cord_risk$threshold, "%"), paste0("≤", cord_risk$threshold, "%")), title = "ASSIGN cardiovascular risk score", legend = "bottom")

fit <- survfit(Surv(p[[times[which_outcome]]], p[[outcomes[which_outcome]]]) ~ cac_med, data = p)
p_2 <- ggsurvplot(fit, pval = TRUE, pval.method = TRUE, risk.table = FALSE, test.for.trend = FALSE, conf.int = FALSE,
                  xlab = "Years", palette = c("#00a1aa", "#aa004c"), ylim = c(0.87,1), pval.coord = c(2, 0.87), pval.method.coord = c(0,0.87),
                  legend.labs = c(paste0(">", round(cord_cac$threshold, 2), "%"), paste0("≤", round(cord_cac$threshold, 2), "%")), title = "Agatston-score", legend = "bottom")

fit <- survfit(Surv(p[[times[which_outcome]]], p[[outcomes[which_outcome]]]) ~ obs_med, data = p)
p_3 <- ggsurvplot(fit, pval = TRUE, pval.method = TRUE, risk.table = FALSE, test.for.trend = FALSE, conf.int = FALSE,
                  xlab = "Years", palette = c("#aa004c", "#00a1aa"), ylim = c(0.87,1), pval.coord = c(2, 0.87), pval.method.coord = c(0,0.87),
                  legend.labs = c("Obstructive CAD", "Non-obstructive CAD"), title = "Obstructive CAD", legend = "bottom")



fit <- survfit(Surv(p[[times[which_outcome]]], p[[outcomes[which_outcome]]]) ~ ncppb_med, data = p)
p_4 <- ggsurvplot(fit, pval = TRUE, pval.method = TRUE, risk.table = FALSE, test.for.trend = FALSE, conf.int = FALSE,
                  xlab = "Years", palette = c("#00a1aa", "#aa004c"), ylim = c(0.87,1), pval.coord = c(2, 0.87), pval.method.coord = c(0,0.87),
                  legend.labs = c(paste0(">", round(cord_ncppb$threshold, 2), "%"), paste0("≤", round(cord_ncppb$threshold, 2), "%")), title = "Non-calcified plaque burden", legend = "bottom")

fit <- survfit(Surv(p[[times[which_outcome]]], p[[outcomes[which_outcome]]]) ~ cppb_med, data = p)
p_5 <- ggsurvplot(fit, pval = TRUE, pval.method = TRUE, risk.table = FALSE, test.for.trend = FALSE, conf.int = FALSE,
                  xlab = "Years", palette = c("#00a1aa", "#aa004c"), ylim = c(0.87,1), pval.coord = c(2, 0.87), pval.method.coord = c(0,0.87),
                  legend.labs = c(paste0(">", round(cord_cppb$threshold, 2), "%"), paste0("≤", round(cord_cppb$threshold, 2), "%")), title = "Calcified plaque burden", legend = "bottom")

fit <- survfit(Surv(p[[times[which_outcome]]], p[[outcomes[which_outcome]]]) ~ lancppb_med, data = p)
p_6 <- ggsurvplot(fit, pval = TRUE, pval.method = TRUE, risk.table = FALSE, test.for.trend = FALSE, conf.int = FALSE,
                  xlab = "Years", palette = c("#00a1aa", "#aa004c"), ylim = c(0.87,1), pval.coord = c(2, 0.87), pval.method.coord = c(0,0.87),
                  legend.labs = c(paste0(">", round(cord_lancppb$threshold, 2), "%"), paste0("≤", round(cord_lancppb$threshold, 2), "%")), title = "Low-attenuation plaque burden", legend = "bottom")



fit <- survfit(Surv(p[[times[which_outcome]]], p[[outcomes[which_outcome]]]) ~ magenta_norm_NCP_med, data = p)
p_7 <- ggsurvplot(fit, pval = TRUE, pval.method = TRUE, risk.table = FALSE, test.for.trend = FALSE, conf.int = FALSE,
                  xlab = "Years", palette = c("#aa004c", "#00a1aa"), ylim = c(0.87,1), pval.coord = c(2, 0.87), pval.method.coord = c(0,0.87),
                  legend.labs = c(paste0(">", round(cord_magenta_norm_NCP$threshold, 2)), paste0("≤", round(cord_magenta_norm_NCP$threshold, 2))), title = "Eigen-NCP 5", legend = "bottom")

fit <- survfit(Surv(p[[times[which_outcome]]], p[[outcomes[which_outcome]]]) ~ brown_norm_CP_med, data = p)
p_8 <- ggsurvplot(fit, pval = TRUE, pval.method = TRUE, risk.table = FALSE, test.for.trend = FALSE, conf.int = FALSE,
                  xlab = "Years", palette = c("#aa004c", "#00a1aa"), ylim = c(0.87,1), pval.coord = c(2, 0.87), pval.method.coord = c(0,0.87),
                  legend.labs = c(paste0(">", round(cord_brown_norm_CP$threshold, 2)), paste0("≤", round(cord_brown_norm_CP$threshold, 2))), title = "Eigen-CP 1", legend = "bottom")

p_clinical <- ggarrange(p_1$plot, p_2$plot, p_3$plot, ncol = 3, nrow = 1)
p_pb  <- ggarrange(p_4$plot, p_5$plot, p_6$plot, ncol = 3, nrow = 1)
p_rad <- ggarrange(p_7$plot, p_8$plot, ncol = 3, nrow = 1)
all_p <- ggarrange(p_clinical, p_pb, p_rad,  ncol = 1, nrow = 3, heights = c(1,1,1))

ggplot2::ggsave(path = folder_img, filename = "KM.pdf", plot = all_p, device = "pdf", width = 12, height = 12, units = "in")