#####
# PROGRAM CODE FOR ARTICLE: SCOTHEART Radiomics data
# (c) Márton Kolossváry, 2022

# DESCRIPTION: Cox analysis
#####

# INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load packages =====
library(data.table); library(survival)

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
prefix <- "COX_Vol_FINAL_"
load(paste0(folder_cache, "FINAL_", which_type, "_", which_d, ".RData"))


## Parameters -----
log_scale <- TRUE  # Should plaque burden covariates be log scaled (TRUE)
limit_CAD <- FALSE # Limit data to only those with CAD (sub-analyses)
rad_bin   <- NULL  # Binarize radiomics to given extreme values (NULL)


covariates <- list(NULL,
                   # c("ncppb", "cppb", "lancppb", "risk", "cac", "obs")
                   c("Total.NCP.vol..mm3.", "Total.CP.vol..mm3.", "Total.LD.NCP.volume..mm3.", "risk", "cac", "obs")
                   )

outcomes <- c("MI_cardiac_death_observed_last")
times    <- c("time_MI_cardiac_death_observed_last")

predictors <- paste0(colnames(p)[grep("_norm", colnames(p),fixed = TRUE)], "")
p$Total.NCP.vol..mm3. <- p$Total.NCP.vol..mm3./100
p$Total.CP.vol..mm3. <- p$Total.CP.vol..mm3./100
p$Total.LD.NCP.volume..mm3. <- p$Total.LD.NCP.volume..mm3./100

plaque_var <- c("pb", "cppb", "ncppb", "lancppb")
# =============================================================================


# RUN COX MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output <- NULL

for(predictor in predictors) {
 for(outcome_i in 1:length(outcomes)) {
  for(covariate in covariates) {
   tryCatch({
    if(is.null(unlist(covariate))) {
     mod <- survival::coxph(
      formula(paste0("survival::Surv(time =", times[outcome_i], ", event=", outcomes[outcome_i], ") ~ ", predictor)), data = p)
    } else {
     mod <- survival::coxph(
      formula(paste0("survival::Surv(time =", times[outcome_i], ", event=", outcomes[outcome_i], ") ~ ", predictor, "+", paste0(unlist(covariate), collapse = "+"))), data = p)
    }}, error = function(e) {})
   
   
   summ   <- broomExtra::tidy(mod, exponentiate = TRUE, conf.int = TRUE)[1, ] # Change if polynominal model
   summ   <- cbind(outcomes[outcome_i], paste0(as.character(unlist(covariate)), collapse = " "), summ)
   summ$n_cases <- dim(p[get(outcomes[outcome_i]) == 1])[1]
   summ$n_pop   <- dim(p)[1]
   
   output <- rbind(output, summ)
  }
 }
}
colnames(output) <- c("Outcome", "Covariates", "EigenRadiomics", "HR", "SE", "Statistic", "p", "ConfLow", "ConfHigh", "nCases", "nPopuliation")
write.csv(output, paste0(folder_stat, prefix, which_type, "_", which_d, ".csv"), row.names = FALSE)

rm(list = ls()); gc(); dev.off()