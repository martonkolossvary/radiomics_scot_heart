#####
# PROGRAM CODE FOR ARTICLE: SCOTHEART Radiomics data
# (c) Márton Kolossváry, 2023

# DESCRIPTION: Table 1

# INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load packages =====
invisible(lapply(list("data.table", "gtsummary", "ggstatsplot", "flextable", "officer"), require, character.only = TRUE))


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
prefix <- "COX_log_"
which_outcome <- 1

load(paste0(folder_cache, "FINAL_ALL_", which_d, ".RData"))
p_eigen <- colnames(p)[grep("_norm", colnames(p),fixed = TRUE)]
p_eigen <- p_eigen[-length(p_eigen)]
eigen_all <- p[, ..p_eigen]
colnames(eigen_all) <- paste0(colnames(eigen_all), "_ALL")
p[p == "NA"] <- NA
p$MI_cardiac_death_observed_last <- as.logical(p$MI_cardiac_death_observed_last)
# SELECT VARIABLES ====================
v <- c("ager", "gender", "bmir", "hypertension", "diabetesmellitusr", "smokinghabit", "chdfamilyhistory", "totalcholestorol", "hdlcholestorol",
       "chestpaindiagnosis", "risk", "cacscore", "cac", "hasCAD", "obs",
       "Total.plaque.vol..mm3.", "Total.NCP.vol..mm3.", "Total.CP.vol..mm3.", "Total.LD.NCP.volume..mm3.",
       "Total.TP.burden", "Total.NCP.burden", "Total.CP.burden", "Total.LD.NCP.burden",
       "pb", "ncppb", "cppb", "lancppb")
v_names <- c("Age [y]", "Sex", "BMI [kg/m2]", "Hypertension", "Diabetes", "Smoking status", "Family history of heart disease", "Total cholesterol [mmol/L]", "HDL cholesterol [mmol/L]",
             "Type of chest pain", "ASSIGN cardiovascular risk [%]", "Agatston-score [units]", "Agatston-score [log2 units]", "Any plaque on CTA", "Obstructive CAD on CTA",
             "Total plaque volume [mm3]", "Noncalcified plaque volume [mm3]", "Calcified plaque volume [mm3]", "Low attenuation noncalcified plaque volume [mm3]",
             "Total-PB [%]", "Noncalcified-PB [%]", "Calcified-PB [%]", "Low attenuation noncalcified-PB [%]",
             "Total-PB [log2]", "Noncalcified-PB [log2]", "Calcified-PB [log2]", "Low attenuation noncalcified-PB [log2]")


# TABLE 1 ====================
## Create PARAMETRIC table -----
t1_all <- tbl_summary(p, by = "MI_cardiac_death_observed_last", statistic = list(all_continuous() ~ "{mean} ± {sd}",
                                                           all_categorical() ~ "{n} ({p}%)"),
                      digits = list(all_continuous() ~ 2),
                      include = any_of(v), sort = list(everything() ~ "alphanumeric"), missing = "no") %>%
 add_overall() %>%
 add_p(list(all_continuous() ~ "t.test", all_categorical() ~ "chisq.test"),
       pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
 modify_caption("Table 1.") %>%
 modify_header(label ~ "Variable")
t1_all$table_body$label[t1_all$table_body$label %in% v] <- v_names
t1_all <- as_flex_table(t1_all)
t1_all <- fontsize(t1_all, size = 8, j = 1:5)
t1_all <- width(t1_all, width = 1.5)
# t1_all <- line_spacing(t1_all, space = 1, part = "all")

save_as_docx(t1_all, path = paste0(folder_stat, "Table_1_parametric.docx"), pr_section = prop_section(page_size = page_size(orient = "landscape")))


