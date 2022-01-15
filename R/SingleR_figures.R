# conda activate r4.1.1
library(Seurat)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame
###############################
pred = readRDS("output/mt_20220114_singleR_azimuth_BlueEncode.rds")
object = readRDS("data/mt_20220114.rds")

singlerDF = data.frame("celltype" = pred$pruned.labels,
                       row.names = rownames(pred))
table(is.na(singlerDF$celltype))
singlerDF$celltype[is.na(singlerDF$celltype)]= "unknown"

