# set defaults & load libraries:
# setwd("C:/Users/vverma3/Desktop/FM_flow_cytometry")
options(device = "RStudioGD", digits = 4, verbose = T)
x <- c(
  "CATALYST", "cowplot", "flowCore", "tidyverse", "diffcyt",
  "SingleCellExperiment", "scater", "readxl", "pheatmap",
  "CytoNorm", "FlowSOM", "gridExtra", "pheatmap", "ggcyto"
)
lapply(x, require, character.only = T) # ~ 19s

# make sce
files <- list.files(
  path = "./dump_neg_fcs",
  pattern = "\\.fcs$", full.names = T
)
PBMC_fs <- read.flowSet(files,
  transformation = F,
  truncate_max_range = F
)
file_name <- as.character(pData(PBMC_fs)$name)

# https://github.com/HelenaLC/CATALYST/issues/103
for (i in 1:length(file_name)) {
  keyword(PBMC_fs@frames[[file_name[i]]])[["$CYT"]] <- "FACS"
}

ID <- gsub("\\.fcs$", "", file_name)
df1 <- data.frame(file_name, ID)
pheno <- read_excel("./FM_pheno.xlsx")
PBMC_md <- merge(df1, pheno, by = "ID", all.x = T)
PBMC_md <- PBMC_md[, c(1:5)]
PBMC_md$date <- NULL
for (i in 1:nrow(PBMC_md)) {
  PBMC_md$date[i] <-
    PBMC_fs@frames[[PBMC_md$file_name[i]]]@description[["$DATE"]]
}
names(PBMC_md) <- c("sample_id", "file_name", "condition", "age", "sex", "date")
PBMC_md$sex <- factor(PBMC_md$sex)
PBMC_md$condition <- factor(PBMC_md$condition,
  levels = c("Control", "Case")
)
PBMC_md <- PBMC_md[, c(2, 1, 3:6)]
PBMC_md$patient_id <- PBMC_md$sample_id
for (i in 1:length(file_name)) {
  cbind(
    print(file_name[i]),
    print(length(PBMC_fs@frames[[file_name[i]]]@parameters@data[["desc"]]))
  )
}
colnames(PBMC_fs) <- gsub("FJComp-", "", gsub("-A$", "", colnames(PBMC_fs)))
fcs_colname <- colnames(PBMC_fs)
antigen <- c(
  "TIGIT", "CD16", "CD57", "CD226", "CD3", "CD56", "CD107a",
  "CD335", "CD159c", "CD158e", "CD314", "CD96", "CD8a", "CD159a"
)
marker_class <- c(
  "state", "type", "type", "state", "type", "type", "state",
  "type", "type", "type", "type", "state", "type", "type"
)
PBMC_panel <- data.frame(fcs_colname, antigen, marker_class,
  stringsAsFactors = F
)
PBMC_panel <- PBMC_panel[c(5, 13, 6, 2, 3, 8, 14, 9, 11, 10, 7, 1, 12, 4), ]
row.names(PBMC_panel) <- NULL
PBMC_panel
sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md,
  md_cols = list(
    file = "file_name", id = "sample_id",
    factors = c("condition", "sex", "age", "date")
  ),
  transform = T, cofactor = 150
)

sce@metadata[["experiment_info"]][["age"]] <- as.numeric(as.character(sce@metadata[["experiment_info"]][["age"]]))
rm(list = (setdiff(ls(), "sce")))
gc()

sce <- cluster(sce,
  features = NULL,
  xdim = 10, ydim = 10, maxK = 50,
  verbose = T, seed = 1
)
gc()
saveRDS(sce, "sce_raw_105.RDS")

# make df for QC: ~
ptm <- proc.time()
ID <- as.character(sce@colData@listData[["sample_id"]])
clust <- as.character(sce@colData@listData[["cluster_id"]])
date <- as.character(sce@colData@listData[["date"]])
df <- data.frame(cbind(ID, clust, date))
meta20 <- sce@metadata[["cluster_codes"]][["meta20"]]
meta50 <- sce@metadata[["cluster_codes"]][["meta50"]]
clust <- c(1:length(meta20))
df1 <- data.frame(clust, meta20, meta50)
df <- merge(df, df1, by = "clust", all.x = T) # ~ 3 min
rm(list = (setdiff(ls(), c("sce", "df"))))
gc()
m <- matrix(
  data = NA, nrow = nrow(df),
  ncol = length(rownames(assay(sce, "exprs"))), byrow = FALSE,
  dimnames = list(1:nrow(df), rownames(assay(sce, "exprs")))
)
for (i in 1:length(rownames(assay(sce, "exprs")))) {
  m[, i] <- assay(sce[rownames(assay(sce, "exprs"))[i], ], "exprs")[1, ]
}
df <- cbind(df, m)
saveRDS(df, "raw_qc_data.RDS")




# rm(list = (setdiff(ls(),"df")));gc()
# sce <- readRDS("sce_raw_105.RDS")
# Do extreme univariate outlier cells belong to same samples?

pheatmap(as.matrix(table(df$CD3, df$ID)),
  color = c("blue", "grey29", "grey", "white", "grey", "grey29", "red"),
  scale = "row",
  angle_col = 90, fontsize = 5, main = "row_scaled"
)





cairo_pdf(filename = "outlier_hm.pdf", width = 15, height = 7, onefile = T)
pheatmap(as.matrix(table(df$clust, df$ID)),
  color = c("blue", "grey29", "grey", "white", "grey", "grey29", "red"),
  scale = "row",
  angle_col = 90, fontsize = 5, main = "row_scaled"
)
pheatmap(as.matrix(table(df$clust, df$ID)),
  color = c("blue", "grey29", "grey", "white", "grey", "grey29", "red"),
  scale = "column",
  angle_col = 90, fontsize = 5, main = "col_scaled"
)
pheatmap(as.matrix(table(df$clust, df$ID)),
  color = c("blue", "grey29", "grey", "white", "grey", "grey29", "red"),
  scale = "none",
  angle_col = 90, fontsize = 5, main = "not_scaled"
)
dev.off()


cairo_pdf(
  filename = "outlier_clust_hm.pdf",
  width = 15, height = 7, onefile = T
)
pheatmap(as.matrix(table(df$clust, df$ID)),
  color = c(
    "blue", "lightblue", "grey70", "grey80",
    "white", "grey80", "grey70", "lightpink1", "red"
  ),
  scale = "row",
  angle_col = 90, fontsize = 5, main = "row_scaled"
)
pheatmap(as.matrix(table(df$clust, df$ID)),
  color = c(
    "blue", "lightblue", "grey70", "grey80",
    "white", "grey80", "grey70", "lightpink1", "red"
  ),
  scale = "column",
  angle_col = 90, fontsize = 5, main = "column_scaled"
)
pheatmap(as.matrix(table(df$ID, df$clust)),
  color = colorRampPalette(c("grey90", "firebrick1"))(5),
  scale = "row",
  angle_col = 90, fontsize = 5, main = "col_scaled"
)
dev.off()

# so far F40 and F107 are forming case-specific clusters
# what about extreme expressions?
p <- plotMedExprs(sce, facet = "antigen")
p$facet$params$ncol <- 4
p



rownames(assay(sce, "exprs"))



CD56 <- assay(sce["CD56", ], "exprs")[1, ]

assay(sce["TIGIT", ], "exprs")[1, ]




clust <- sce@colData@listData[["cluster_id"]]
date <- sce@colData@listData[["date"]]
ID <- sce@metadata[["experiment_info"]][["sample_id"]]
meta20 <- sce@metadata[["cluster_codes"]][["meta20"]]
clust <- c(1:length(meta20))
df1 <- data.frame(clust, meta20)
df2 <- merge(df, df1, by = "clust", all.x = T)
chisq.test(table(df2$date, df2$meta20))
# X-squared = 928266, df = 57, p-value <2e-16
m <- as.matrix(table(df$meta20, df$date))
cairo_pdf(filename = "clust_batch_hm.pdf", width = 5, height = 5)
pheatmap(m) # outlier: clust 3 / 28-AUG batch
dev.off()








# marker-specific outliers:

cairo_pdf(filename = "outlier_marker.pdf", width = 8, height = 15, onefile = T)
plotExprHeatmap(sce,
  row_anno = F, scale = "last", row_clust = T, col_clust = F,
  row_dend = T, col_dend = T, bin_anno = F, bars = F, perc = F,
  hm_pal = c("blue", "grey", "white", "grey", "red")
)
dev.off()
# E28, F107 - extreme expressions

cairo_pdf(filename = "outlier_hm.pdf", width = 10, height = 7, onefile = T)
pheatmap(as.matrix(table(cluster_ids(sce, "meta20"), sample_ids(sce))),
  color = c("blue", "grey", "white", "grey", "red"), scale = "row",
  angle_col = 90, main = "rows scaled"
)
pheatmap(as.matrix(table(cluster_ids(sce, "meta5"), sample_ids(sce))),
  color = c("blue", "grey", "white", "grey", "red"), scale = "column",
  angle_col = 90, main = "columns scaled"
)
dev.off()

# column: F107, E28, E58, F45
# row: F107, E58, E40, F96, F40

# E28 (cntrl) and F107 (case) were dropped:
dir.create("./dump_neg_fcs/bad_fcs")
file.copy(
  c("./dump_neg_fcs/E28.fcs", "./dump_neg_fcs/F107.fcs"),
  "./dump_neg_fcs/bad_fcs/"
)
file.remove(c("./dump_neg_fcs/E28.fcs", "./dump_neg_fcs/F107.fcs"))

# is there a need for batch normalization?
# ------------------------------------------------------------ 0 min
files <- list.files(
  path = "./dump_neg_fcs",
  pattern = "\\.fcs$", full.names = T
)
PBMC_fs <- read.flowSet(files,
  transformation = F,
  truncate_max_range = F
)
file_name <- as.character(pData(PBMC_fs)$name)
ID <- gsub("\\.fcs$", "", file_name)
df1 <- data.frame(file_name, ID)
pheno <- read_excel("./FM_pheno.xlsx")
PBMC_md <- merge(df1, pheno, by = "ID", all.x = T)
PBMC_md <- PBMC_md[, c(1:5)]
PBMC_md$date <- NULL
for (i in 1:nrow(PBMC_md)) {
  PBMC_md$date[i] <-
    PBMC_fs@frames[[PBMC_md$file_name[i]]]@description[["$DATE"]]
}
names(PBMC_md) <- c("sample_id", "file_name", "condition", "age", "sex", "date")
PBMC_md$sex <- factor(PBMC_md$sex)
PBMC_md$condition <- factor(PBMC_md$condition,
  levels = c("Control", "Case")
)
PBMC_md <- PBMC_md[, c(2, 1, 3:6)]
PBMC_md$patient_id <- PBMC_md$sample_id
colnames(PBMC_fs) <- gsub("FJComp-", "", gsub("-A$", "", colnames(PBMC_fs)))
fcs_colname <- colnames(PBMC_fs)
antigen <- c(
  "TIGIT", "CD16", "CD57", "CD226", "CD3", "CD56", "CD107a",
  "CD335", "CD159c", "CD158e", "CD314", "CD96", "CD8a", "CD159a"
)
marker_class <- c(
  "state", "type", "type", "state", "type", "type", "state",
  "type", "type", "type", "type", "state", "type", "type"
)
PBMC_panel <- data.frame(fcs_colname, antigen, marker_class,
  stringsAsFactors = F
)
PBMC_panel
sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md,
  md_cols = list(
    file = "file_name", id = "sample_id",
    factors = c("condition", "sex", "age", "date")
  ),
  transform = T, cofactor = 150
)

sce@metadata[["experiment_info"]][["age"]] <- as.numeric(as.character(sce@metadata[["experiment_info"]][["age"]]))
rm(list = (setdiff(ls(), "sce")))
gc()

sce <- cluster(sce,
  features = NULL,
  xdim = 10, ydim = 10, maxK = 20,
  verbose = T, seed = 1
)
gc()
# ------------------------------------------------------------ 11 min
date <- as.factor(sce@metadata[["experiment_info"]][["date"]])
fm <- as.factor(sce@metadata[["experiment_info"]][["condition"]])
chisq.test(table(date, fm))
# X-squared = 1.9, df = 3, p-value = 0.6
clust <- sce@colData@listData[["cluster_id"]]
date <- sce@colData@listData[["date"]]
chisq.test(table(date, clust))
# X-squared = 1744712, df = 297, p-value <2e-16 (for 100 clust)
df <- data.frame(clust, date)
meta20 <- sce@metadata[["cluster_codes"]][["meta20"]]
clust <- c(1:length(meta20))
df1 <- data.frame(clust, meta20)
df2 <- merge(df, df1, by = "clust", all.x = T)
chisq.test(table(df2$date, df2$meta20))
# X-squared = 928266, df = 57, p-value <2e-16
m <- as.matrix(table(df2$meta20, df2$date))
cairo_pdf(filename = "clust_batch_hm.pdf", width = 5, height = 5)
pheatmap(m) # outlier: clust 3 / 28-AUG batch
dev.off()

saveRDS(sce, "sce_raw_103.RDS")

system.time(sce <- readRDS("sce_raw_103.RDS")) # ~ 70 sec

cairo_pdf(filename = "outlier_hm1.pdf", width = 10, height = 7, onefile = T)
pheatmap(as.matrix(table(cluster_ids(sce, "meta20"), sample_ids(sce))),
  color = c("blue", "grey", "white", "grey", "red"), scale = "row",
  angle_col = 90, main = "rows scaled"
)
pheatmap(as.matrix(table(cluster_ids(sce, "meta20"), sample_ids(sce))),
  color = c("blue", "grey", "white", "grey", "red"), scale = "column",
  angle_col = 90, main = "columns scaled"
)
dev.off()



### *************** Evidence of batch-effect found !! *************** ###
##### Cytonorm #####
rm(list = ls())
gc()
data <- read.delim("S:/FM-NK-Flowcytometry/batch_data.txt",
  stringsAsFactors = F
)
data$Type <- as.character(data$Type)
data$Batch <- as.factor(data$Batch)
data$Type <- c("1" = "Train", "2" = "Validation")[data$Type]
table(data$Type, data$Batch)
#            A  B  C  D
# Train       1  1  1  1
# Validation 20 17 49 13

train_data <- dplyr::filter(data, Type == "Train")
validation_data <- dplyr::filter(data, Type == "Validation")
ff <- flowCore::read.FCS(data$Path[1])
channels <- flowCore::colnames(ff)
fcTransform <- flowCore::arcsinhTransform(
  transformationId = "fcTransform",
  a = 0, b = (1 / 150), c = 0
)
fcTransform.reverse <- function(x) {
  return(sinh(x) / (1 / 150))
}

transformList <- flowCore::transformList(channels, fcTransform)
transformList.reverse <- flowCore::transformList(
  channels,
  fcTransform.reverse
)

fsom <- prepareFlowSOM(train_data$Path, channels,
  nCells = 150000,
  FlowSOM.params = list(
    xdim = 20, ydim = 20,
    nClus = 50, scale = F
  ),
  transformList = transformList, seed = 1
) # ~ 30 sec
# calculate coeff of variation to find optimum no. of clusters
cvs <- testCV(fsom, cluster_values = 3:50, plot = T, verbose = T) # ~ 15 min

saveRDS(fsom, "cytonorm_fsom.RDS")
saveRDS(cvs, "cytonorm_cvs.RDS")
cvs[["cvs"]][["7"]]
cvs[["cvs"]][["6"]] # all CVs < 1

cvs$pctgs$`25` # percentages

model <- CytoNorm.train(
  files = train_data$Path,
  labels = train_data$Batch,
  channels = channels,
  transformList = transformList,
  FlowSOM.params = list(
    nCells = 150000,
    xdim = 10,
    ydim = 10,
    nClus = 6,
    scale = FALSE
  ),
  normMethod.train = QuantileNorm.train,
  normParams = list(nQ = 101, goal = "mean"),
  seed = 1, verbose = T, plot = T
)
gc()
CytoNorm.normalize(
  model = model,
  files = validation_data$Path,
  labels = validation_data$Batch,
  transformList = transformList,
  transformList.reverse = transformList.reverse,
  normMethod.normalize = QuantileNorm.normalize,
  outputDir = "normalized", prefix = "nrm_",
  clean = T, verbose = T
)






# normalization introduced infinite values for some cells
# we filter those out:
sce <- filterSCE(sce, colSums(is.infinite(assay(sce, "exprs"))) == 0)
# this step takes ~14 min
sce <- cluster(sce,
  features = "type",
  xdim = 10, ydim = 10, maxK = 10,
  verbose = T, seed = 1
)
gc()

saveRDS(sce, "sce_nrm_101.RDS")
sce <- readRDS("sce_nrm_101.RDS")

date <- as.factor(sce@metadata[["experiment_info"]][["date"]])
fm <- as.factor(sce@metadata[["experiment_info"]][["condition"]])
chisq.test(table(date, fm))
# X-squared = 1.7165, df = 3, p-value = 0.6333

clust <- sce@colData@listData[["cluster_id"]]
date <- sce@colData@listData[["date"]]
chisq.test(table(date, clust))
# X-squared = 2422746, df = 297, p-value < 2.2e-16 (for 100 clust)

df <- data.frame(clust, date)
meta10 <- sce@metadata[["cluster_codes"]][["meta10"]]
clust <- c(1:length(meta10))
df1 <- data.frame(clust, meta10)
df2 <- merge(df, df1, by = "clust", all.x = T)
chisq.test(table(df2$date, df2$meta10))
# X-squared = 1148191, df = 27, p-value < 2.2e-16 (chi-sq improved by 219149)
table(df2$meta10, df2$date)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> $$$$ Without F107 $$$$
rm(list = ls())
gc()
files <- list.files(
  path = "./normalized",
  pattern = "\\.fcs$", full.names = T
)
PBMC_fs <- read.flowSet(files,
  transformation = F,
  truncate_max_range = F
)
rm(list = (setdiff(ls(), "PBMC_fs")))
gc()

file_name <- as.character(pData(PBMC_fs)$name)
ID <- gsub("^nrm_", "", gsub("\\.fcs$", "", file_name))
df1 <- data.frame(file_name, ID)
pheno <- read_excel("./FM_pheno.xlsx")
PBMC_md <- merge(df1, pheno, by = "ID", all.x = T)
PBMC_md <- PBMC_md[, c(1:5)]
PBMC_md$date <- NULL
for (i in 1:nrow(PBMC_md)) {
  PBMC_md$date[i] <-
    PBMC_fs@frames[[PBMC_md$file_name[i]]]@description[["$DATE"]]
}
names(PBMC_md) <- c("sample_id", "file_name", "condition", "age", "sex", "date")
PBMC_md$sex <- factor(PBMC_md$sex)
PBMC_md$condition <- factor(PBMC_md$condition,
  levels = c("Control", "Case")
)
PBMC_md <- PBMC_md[, c(2, 1, 3:6)]
PBMC_md$patient_id <- PBMC_md$sample_id
rm(list = (setdiff(ls(), c("PBMC_fs", "PBMC_md"))))
gc()
table(PBMC_md$condition)
# Control    Case
#      54      46
colnames(PBMC_fs) <- gsub("FJComp-", "", gsub("-A$", "", colnames(PBMC_fs)))
fcs_colname <- colnames(PBMC_fs)
antigen <- c(
  "TIGIT", "CD16", "CD57", "CD226", "CD3", "CD56", "CD107a",
  "CD335", "CD159c", "CD158e", "CD314", "CD96", "CD8a", "CD159a"
)
marker_class <- c(
  "state", "type", "type", "state", "type", "type", "state",
  "type", "type", "type", "type", "state", "type", "type"
)
PBMC_panel <- data.frame(fcs_colname, antigen, marker_class,
  stringsAsFactors = F
)
PBMC_panel
rm(list = (setdiff(ls(), c("PBMC_fs", "PBMC_md", "PBMC_panel"))))
gc()

# this step takes ~1 min (local)
sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md,
  md_cols = list(
    file = "file_name", id = "sample_id",
    factors = c("condition", "sex", "age", "date")
  ),
  transform = T, cofactor = 150
)
# normalization introduced infinite values for some cells:
table(colSums(is.infinite(assay(sce, "exprs"))))
#        0        1
# 14999921       79
# we filter those out:
sce <- filterSCE(sce, colSums(is.infinite(assay(sce, "exprs"))) == 0)

rm(list = (setdiff(ls(), "sce")))
gc()

# this step takes ~14 min

sce <- cluster(sce,
  features = "type",
  xdim = 10, ydim = 10, maxK = 20,
  verbose = T, seed = 1
)
gc()
# fix age from factor to numeric:
sce@metadata[["experiment_info"]][["age"]] <- as.numeric(sce@metadata[["experiment_info"]][["age"]])
saveRDS(sce, "sce_nrm_100.RDS")
# sce <- readRDS("sce_nrm_100.RDS")

date <- as.factor(sce@metadata[["experiment_info"]][["date"]])
fm <- as.factor(sce@metadata[["experiment_info"]][["condition"]])
chisq.test(table(date, fm))
# X-squared = 2.1489, df = 3, p-value = 0.5421

clust <- sce@colData@listData[["cluster_id"]]
date <- sce@colData@listData[["date"]]
chisq.test(table(date, clust))
# X-squared = 1926493, df = 297, p-value < 2.2e-16

df <- data.frame(clust, date)
meta10 <- sce@metadata[["cluster_codes"]][["meta10"]]
clust <- c(1:length(meta10))
df1 <- data.frame(clust, meta10)
df2 <- merge(df, df1, by = "clust", all.x = T)
chisq.test(table(df2$date, df2$meta10))
# X-squared = 615097, df = 27, p-value < 2.2e-16
#  (chi-sq further reduced by 533094)
table(df2$meta10, df2$date)

design <- createDesignMatrix(ei(sce),
  cols_design = c("condition", "age", "sex", "date")
)
contrast <- createContrast(c(0, 1, 0, 0, 0, 0, 0))
nrow(contrast) == ncol(design) # TRUE
data.frame(parameters = colnames(design), contrast)


### Step 6: diffcyt
res_DA <- diffcyt(sce,
  clustering_to_use = "meta20",
  analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
  design = design, contrast = contrast, verbose = T
)
topTable(res_DA, format_vals = T)

plotAbundances(sce,
  k = "meta20", by = "cluster_id",
  group_by = "condition"
)

plotAbundances(sce, k = "meta20", by = "sample_id", group_by = "condition")

# subset type-marker expression matrix
es <- assay(sce, "exprs")

# run diffcyt with formula:

ei <- metadata(sce)$experiment_info
ei$age <- as.numeric(sce@metadata[["experiment_info"]][["age"]])
da_formula1 <- createFormula(ei,
  cols_fixed = c("condition", "age", "sex"),
  cols_random = "date"
)
contrast <- createContrast(c(0, 1, 0, 0))
da_res1 <- diffcyt(sce,
  formula = da_formula1, contrast = contrast,
  analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
  clustering_to_use = "meta10", verbose = T
)
topTable(da_res1, format_vals = T)



form <- createFormula(PBMC_md, cols_fixed = c("condition", "age", "sex", "date"))

res_DA1 <- diffcyt(sce,
  clustering_to_use = "meta20",
  analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
  formula = form, design = design,
  contrast = contrast, verbose = T
)
topTable(res_DA1, format_vals = T)

# what about other outliers:
sce <- readRDS("sce_nrm_100.RDS")

plotAbundances(sce,
  k = "meta10", by = "cluster_id",
  group_by = "condition"
)



x <- as.matrix(table(cluster_ids(sce, "meta10"), sample_ids(sce)))
pheatmap(as.matrix(table(cluster_ids(sce, "meta10"), sample_ids(sce))), scale = "row", fontsize = 8)

plotClusterExprs(sce, k = "meta10", features = "CD8a")

es <- assay(sce["CD56", ], "exprs")
es <- as.numeric(es[1, ])


plotExprHeatmap(sce,
  features = "CD8a", by = "both", k = "som100",
  scale = "never", col_clust = FALSE, row_anno = FALSE, bars = FALSE
)
