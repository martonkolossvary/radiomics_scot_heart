#####
# PROGRAM CODE FOR ARTICLE: SCOTHEART Radiomics data
# (c) Márton Kolossváry, 2023

# DESCRIPTION: Run WGCNA and create eigen radiomic features
#####


# INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load packages =====
library(data.table); library(WGCNA)

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
first_rad  <- 8 #First radiomic parameter column number

### Remove those who have plaque but no nifti
p <- readxl::read_xlsx(paste0(folder_data, "SCOT-HEART Per-Patient Masterlist.xlsx"))
setDT(p)
load(paste0(folder_cache, "DATA_", "ALL", "_", which_d, "_distmat.RData"))
exc <- p[(!subjectno %in% d$patient_ID) & Total.plaque.vol..mm3. >0]$subjectno
p   <- p[!subjectno %in% exc]

load(paste0(folder_cache, "DATA_", which_type, "_", which_d, "_distmat.RData"))
load(paste0(folder_cache, "dist_mat_", which_type, "_", which_d, ".RData"))


## Parameters -----
POWERS   <- 2:10 # Powers used for scale free property analysis
COR_FUNC <- "pearson" # Correlation type cor (corOptions sets p) or bicor if similarity was not calculated previously
BETA     <- 7  # ALL_EP:7 NCP_EP:5  CP_EP:8
nSets    <- 1 # Number of temporal sets

plotCols    <- c(2,5,6,7,9,10) # Which statistics to plot
colNames    <- c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
                 "Max connectivity", "Network Density", "Network Heterogeneity") # Names of statistics
colors      <- ggsci::pal_npg()(2)
d_names     <- c("Pearson correlation", "Biweight midcorrelation")
# =============================================================================


# Scale-free plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Do scale free property analysis using precalculated similarity matrix -----
powerTables <- vector(mode = "list", length = 1); names(powerTables) <- "Similarity"
powerTables <- lapply(1, function(x) {
 pickSoftThreshold.fromSimilarity(similarity = m_dist, powerVector = POWERS, verbose = 2, moreNetworkConcepts = TRUE, removeFirst = TRUE, blockSize = NULL)})

## Get plot y ranges -----
ylim <-  matrix(NA, nrow = 2, ncol = length(plotCols))
for (set in 1) {
 for (col in 1:length(plotCols)) {
  ylim[1, col] = min(ylim[1, col], powerTables[[set]]$fitIndices[, plotCols[col]], na.rm = TRUE)
  ylim[2, col] = max(ylim[2, col], powerTables[[set]]$fitIndices[, plotCols[col]], na.rm = TRUE)
 }
}

## Plot network statistics -----
pdf(file = paste0(folder_img, "Network_statistics_", which_type, "_", which_d, ".pdf"), width = 12, height = 6)
par(mfcol = c(2,3)); par(mar = c(4.2, 4.2 , 2.2, 4.5), xpd = TRUE); cex1 = 0.7

for (col in 1:(length(plotCols))) {
 for (set in 1) {
  if (set==1) {
   plot(powerTables[[set]]$fitIndices[,1], (1*powerTables[[set]]$fitIndices[,3]<0)*powerTables[[set]]$fitIndices[,2],
        xlab="Soft Threshold (power)", ylab = colNames[col], type = "n", ylim = ylim[, col], main = colNames[col], xaxt = "n")
   axis(side = 1, at=c(seq(min(POWERS), max(POWERS), 2))) # Set X axis
   addGrid() }
  if (col==1) {
   lines(powerTables[[set]]$fitIndices[,1], (1*powerTables[[set]]$fitIndices[,3]<0)*powerTables[[set]]$fitIndices[,2],
         type = "l", col = colors[set])
  } else {
   lines(powerTables[[set]]$fitIndices[,1], powerTables[[set]]$fitIndices[, plotCols[col]],
         type = "l", col = colors[set])
  }
  legend("topright", legend = "Similarity", col = colors, pch = 20, inset=c(0,0))
 }
}
dev.off()
# =============================================================================


# Calculate WGCNA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Create radiomics data table for individuals without plaque -----
p0 <- rbindlist(list(d[0, ], data.table(patient_ID = as.character(p[!(subjectno %in% d$patient_ID)]$subjectno))), fill = TRUE)
p0[, first_rad:(dim(d)[2])][is.na(p0[, first_rad:(dim(d)[2])])] <- 0
p0$plaque_ID <- "NONE"
p0$Group_name <- paste0(p0$patient_ID, "__", p0$plaque_ID, ".nii")
d_WGCNA <- rbindlist(list(d, p0))

multiExpr <- list(); multiExpr$BAS$data <- d_WGCNA[, first_rad:(dim(d)[2])] # change between d / d_WGCNA for only plaque / including 0s
nGenes <- dim(m_dist)[1]

## Create consensus TOM ----
adjacencies <- array(0, dim = c(nSets, nGenes, nGenes))
TOM         <- array(0, dim = c(nSets, nGenes, nGenes))
for (set in 1:nSets) {
 adjacencies[set, , ] <- m_dist^BETA
 TOM[set, , ] <- TOMsimilarity(adjacencies[set, , ])
 #Calibration and consensus TOM creation if nSets>1
}
consensusTOM <- pmin(TOM[1, , ])

## Create clustering -----
minModuleSize <- 20

consTree       <- fastcluster::hclust(as.dist(1-consensusTOM), method = "average")
unmergedLabels <- cutreeDynamic(dendro = consTree, distM = 1-consensusTOM, method = "hybrid",
                                deepSplit = 2, cutHeight = 0.995, minClusterSize = minModuleSize,
                                assumeSimpleExternalSpecification = FALSE,
                                pamStage = TRUE, pamRespectsDendro = TRUE) #Same defaults as blockwiseConsensusModules and cutreeDynamic
unmergedColors <- labels2colors(unmergedLabels)
unmergedMEs    <- multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors, excludeGrey = TRUE)

### Merge clusters using only eigengene similarity -----
merge <- mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.10, verbose = 3, corFnc = bicor)
moduleLabels <- merge$newMEs$BAS$validColors
moduleColors <- labels2colors(moduleLabels)
consMEs <- merge$newMEs
# =============================================================================
write.csv(table(moduleColors, dict_rad$rad_class), "/Users/mjk2/Downloads/NCP.csv")

# PLOT WGCNA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Plot dendrogram =====
pdf(file = paste0(folder_img, "WGCNA_dendrogram_", which_type, "_", which_d, ".pdf"), width = 12, height = 12)
plotDendroAndColors(dendro = consTree, 
                    colors = cbind(moduleColors),
                    groupLabels = c("Eigen modules"),
                    cex.dendroLabels = 0.4,
                    dendroLabels = FALSE, hang = 0.05,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus eigen module assignment of radiomic features")
dev.off()

## Plot heatmap =====
pdf(file = paste0(folder_img, "WGCNA_TOM_Heatmap_", which_type, "_", which_d, ".pdf"), width = 12, height = 12)
p_TOMcon <- TOMplot(dissim = consensusTOM, dendro = consTree, Colors = moduleColors, terrainColors = FALSE)
dev.off()

## Plot MDS plot =====
library(plotly)

mds <- cmdscale(1-consensusTOM, k = 3)
module_colors <- factor(moduleColors, levels = unique(moduleColors)) # con_net_all$colors or con_net_all_sp3$color
mds <- cbind(as.data.frame(mds), module_colors); colnames(mds) <- c("x", "y", "z", "color")

## Plot parameters =====
l <- list(font = list(family = "helvetica",size = 14),
          title = list(text='<b> Radiomic Eigen Modules </b>'),
          itemsizing = list(values = "constant",itemwidth = 50),
          orientation = 'h')

s <- list(camera = list(eye = list(x = 0.5, y = 1.5, z = 0.8)), #-1.5
          xaxis = list(title = "Scaling Dimension 1",
                       font = list(family = "helvetica", size = 12),
                       gridcolor = toRGB("black"), gridwidth = 4),
          yaxis = list(title = "Scaling Dimension 2",
                       font = list(family = "helvetica", size = 12),
                       gridcolor = toRGB("black"), gridwidth = 4),
          zaxis = list(title = "Scaling Dimension 3",
                       font = list(family = "helvetica", size = 12),
                       gridcolor = toRGB("black"), gridwidth = 4))
m <- list(l = 0, r = 0, b = 0, t = 0, pad = 1)

fig <- plotly::plot_ly(type = "scatter3d", x = mds$x, y = mds$y, z = mds$z,
                       marker = list(size = 8, opacity = 1),
                       color = mds$color, colors = gplots::col2hex(unique(mds$color)))
fig <- fig %>% layout(legend = NULL, scene = s, margin = m)

setwd(folder_img)
orca(fig, paste0("MDS_plot_", which_type, "_", which_d, ".png"), format = "png", width = 6*300, height = 6*300)
setwd(folder_wd)


# FORMAT AND SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Create final plaque based data -----
d_eigen <- consMEs$BAS$data # data = Module eigengenes; averageExpr = average normalized expression
colnames(d_eigen) <- labels2colors(as.numeric(gsub("ME", "", colnames(d_eigen))))
d_eigen <- cbind(d_WGCNA[, c("patient_ID", "plaque_ID")], d_eigen)

d_final <- merge(d_eigen, d_WGCNA, by = c("patient_ID", "plaque_ID"), all = TRUE)

save(list = c("d_final"),
     file = paste0(folder_cache, "Plaque_", which_type, "_", which_d, ".RData"), envir = .GlobalEnv)
rm(list = ls()); gc(); dev.off()
## ============================================================================