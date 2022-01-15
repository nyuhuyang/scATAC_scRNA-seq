#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
invisible(lapply(c("Seurat","SingleR","SingleCellExperiment",
                   "magrittr","data.table","Matrix"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# Need 64GB ?
set.seed(101)


# ====== load single cell =============
object = readRDS("data/mt_20220114.rds")
sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# ====== load reference =============
blue_encode <- BlueprintEncodeData()
rownames(blue_encode) %<>% tolower %>% Hmisc::capitalize()
# ====== conbime data =============

common <- Reduce(intersect, list(rownames(sce),
                                 rownames(blue_encode)
))
length(common)
table(blue_encode$label.fine)
system.time(trained <- trainSingleR(ref = blue_encode[common,],
                                    labels=blue_encode$label.fine))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = paste0("output/mt_20220114_singleR_azimuth_BlueEncode.rds"))
