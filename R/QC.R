########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
# conda activate r4.1.1
invisible(lapply(c("Seurat","dplyr","ggplot2","scater","magrittr","pbapply","scales",
                   "cowplot"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

if(!dir.exists("data")) dir.create("data")
if(!dir.exists("doc")) dir.create("doc")
if(!dir.exists("output")) dir.create("output")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20210114_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples %<>% filter(seq == "scRNA")
# check missing data
current <- list.files("data",recursive = TRUE) %>%
                gsub("/outs.*","",.) %>%
                gsub(".*/","",.) %>% unique
(current <- current[!grepl(".Rda|RData",current)])
print(paste("missing data:",missing_data <- df_samples$sample.id[!(df_samples$sample.id %in% current)]))
df_samples 
#======1.1.2 record data quality before removing low quanlity cells =========================
# seq == "scRNA"
message("read metrics_summary")
QC_list <- lapply(df_samples$sample.id[df_samples$seq == "scRNA"], function(x){
        tmp = read.csv(file = paste0("data/Processed_scRNA_wTDT_2/",x,
                               "/outs/metrics_summary.csv"))
        t(tmp)
})
names(QC_list) = df_samples$sample.id[df_samples$seq == "scRNA"]

cbind.fill <- function(...){
        nm <- list(...) 
        nm <- lapply(nm, as.matrix)
        n <- max(sapply(nm, nrow)) 
        do.call(cbind, lapply(nm, function (x) 
                rbind(x, matrix(NA, n-nrow(x), ncol(x))))) 
}
QC <- do.call(cbind.fill, QC_list)
colnames(QC) = df_samples$sample.id[df_samples$seq == "scRNA"]
QC[is.na(QC)] = ""

QC["Estimated.Number.of.Cells",] %>% gsub(",","",.) %>% as.numeric %>% sum

write.csv(QC,paste0(path,"metrics_summary.csv"))
df_samples %<>% cbind(t(QC))
rownames(df_samples) = df_samples$sample
openxlsx::write.xlsx(df_samples, file =  paste0(path,"20210114_scRNAseq_info.xlsx"),
                     colNames = TRUE,rowNames = T,borders = "surrounding",
                     colWidths = c(NA, "auto", "auto"), overwrite = TRUE)

## Load the GEX dataset
message("Loading the datasets")
Seurat_list <- list()
for(i in seq_along(df_samples$sample.id)){
        print(s <- df_samples$sample.id[i])
        tmp <- Read10X_h5(filename = paste0("data/Processed_scRNA_wTDT_2/",as.character(s),"/outs/filtered_feature_bc_matrix.h5"))
        if(class(tmp) == "list") tmp = tmp$`Gene Expression`
        colnames(tmp) %<>% gsub("-[0-9+]","",.)
        colnames(tmp) %<>% paste0(df_samples$sample[i],"-",.)
        Seurat_list[[s]] = CreateSeuratObject(tmp,min.cells = 0,names.delim = "-",min.features = 0)
        Progress(which(df_samples$sample.id %in% s), length(df_samples$sample.id))
        }
rm(tmp);GC()

#========1.1.3 g1 QC plots before filteration=================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), Seurat_list)
remove(Seurat_list);GC()
object$orig.ident %<>% factor(levels = df_samples$sample)
table(object$orig.ident)

# read and select mitochondial genes
mito <-"^mt-"
message("mito.genes:")

(mito.features <- grep(pattern = mito, x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = mito)
Idents(object) <- "orig.ident"
Idents(object) %<>% factor(levels = df_samples$sample)
g1 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 1, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=8,angle = 90,hjust = 1,vjust = 0.5),
                      legend.position="none",plot.title = element_text(hjust = 0.5))
})
save(g1,file= paste0(path,"g1","_",length(df_samples$sample.id),"_",gsub("-","",Sys.Date()),".Rda"))

#============1.2 scatter ======================
meta_data = object@meta.data
meta_data$discard = FALSE
for(i in seq_along(df_samples$sample)){
        cell = rownames(meta_data)[meta_data$orig.ident %in% df_samples$sample[i]]
        high.mito <- isOutlier(meta_data[cell,"percent.mt"], nmads=3, type="higher")
        low.lib <- isOutlier(log10(meta_data[cell,"nCount_RNA"]), type="lower", nmad=3)
        low.genes <- isOutlier(log10(meta_data[cell,"nFeature_RNA"]), type="lower", nmad=3)
        discard <- high.mito | low.lib | low.genes
        print(data.frame(HighMito= sum(high.mito),LowLib=sum(low.lib),
                         LowNgenes=sum(low.genes),Discard=sum(discard)))
        meta_data[cell,"discard"] = discard
}
object@meta.data = meta_data
meta_data = object@meta.data
meta_data$discard = FALSE
for(i in seq_along(df_samples$sample)){
        cell = rownames(meta_data)[meta_data$orig.ident %in% df_samples$sample[i]]
        high.mito <- isOutlier(meta_data[cell,"percent.mt"], nmads=3, type="higher")
        low.lib <- isOutlier(log10(meta_data[cell,"nCount_RNA"]), type="lower", nmad=3)
        low.genes <- isOutlier(log10(meta_data[cell,"nFeature_RNA"]), type="lower", nmad=3)
        discard <- high.mito | low.lib | low.genes
        print(data.frame(HighMito= sum(high.mito),LowLib=sum(low.lib),
                         LowNgenes=sum(low.genes),Discard=sum(discard)))
        meta_data[cell,"discard"] = discard
}
object@meta.data = meta_data
table(object$orig.ident, object$discard)

object %<>% subset(subset = discard == FALSE
                   &  nFeature_RNA > 200
                   & nCount_RNA > 500
                   & percent.mt < 10
                   )
# FilterCellsgenerate Vlnplot before and after filteration
Idents(object) = "orig.ident"

g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 1, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=8,angle = 90,hjust = 1,vjust = 0.5),
                      legend.position="none",plot.title = element_text(hjust = 0.5))
})
save(g2,file= paste0(path,"g2","_",length(df_samples$sample.id),"_",gsub("-","",Sys.Date()),".Rda"))

jpeg(paste0(path,"S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g1[[1]]+ggtitle("nFeature_RNA before filteration")+
                        scale_y_log10(limits = c(100,10000)),
        g2[[1]]+ggtitle("nFeature_RNA after filteration")+
                        scale_y_log10(limits = c(100,10000)))
dev.off()
jpeg(paste0(path,"S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g1[[2]]+ggtitle("nCount_RNA before filteration")+
                        scale_y_log10(limits = c(500,100000)),
                g2[[2]]+ggtitle("nCount_RNA after filteration")+
                        scale_y_log10(limits = c(500,100000)))
dev.off()
jpeg(paste0(path,"S1_mito.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g1[[3]]+ggtitle("mito % before filteration")+
                        ylim(c(0,50)),
                g2[[3]]+ggtitle("mito % after filteration")+
                        ylim(c(0,50)))
dev.off()

#====
format(object.size(object),unit = "GB")
saveRDS(object, file = "data/mt_20220114.rds")