########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.1.1
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","harmony","sctransform"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


########################################################################
#
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("output/20220114/20210114_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples$mouse %<>% as.character()
#======1.2 load  Seurat =========================
object = readRDS(file = "data/mt_20220114.rds")
meta.data = object@meta.data
for(i in seq_along(df_samples$sample)){
    cells <- meta.data$orig.ident %in% df_samples$sample[i]
    print(df_samples$sample[i])
    print(table(cells))
    meta.data[cells,"sample"] = df_samples$sample[i]
    meta.data[cells,"mouse"] = df_samples$mouse[i]
    meta.data[cells,"group"] = df_samples$group[i]
    meta.data[cells,"estimated.number.of.cells"] = df_samples$estimated.number.of.cells[i]
    meta.data[cells,"mean.reads.per.cell"] = df_samples$mean.reads.per.cell[i]
    meta.data[cells,"number.of.reads"] = df_samples$number.of.reads[i]
    meta.data[cells,"sequencing.saturation"] = df_samples$sequencing.saturation[i]
}
meta.data$sample %<>% factor(levels = df_samples$sample)

table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data
#======1.6 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()

object_list %<>% pblapply(SCTransform,method = "glmGamPoi",vars.to.regress = "percent.mt")
features <- SelectIntegrationFeatures(object.list = object_list)

options(future.globals.maxSize= object.size(object_list)*1.5)
object_list %<>% pblapply(FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
    x
})
anchors <- FindIntegrationAnchors(object.list = object_list,reference = c(1, 2), reduction = "rpca", 
                                  dims = 1:50)
remove(object_list);GC()
# this command creates an 'integrated' data assay
object <- IntegrateData(anchorset = anchors,normalization.method = "SCT", dims = 1:50)
remove(anchors);GC()
format(object.size(object),unit = "GB")
saveRDS(object, file = "data/mt_20220114.rds")

# Determine the ‘dimensionality’ of the dataset  =========
npcs = 100

DefaultAssay(object) <- "RNA"
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = 100, verbose = FALSE)
jpeg(paste0(path,"ElbowPlot_RNA.jpeg"), units="in", width=10, height=7,res=600)
print(ElbowPlot(object,ndims = 100))
dev.off()
object <- JackStraw(object, num.replicate = 20,dims = npcs)
object <- ScoreJackStraw(object, dims = 1:npcs)

for(i in 0:9){
    a = i*10+1; b = (i+1)*10
    jpeg(paste0(path,"JackStrawPlot_",a,"_",b,".jpeg"), units="in", width=10, height=7,res=600)
    print(JackStrawPlot(object, dims = a:b))
    dev.off()
    Progress(i, 9)
}
p.values = object[["pca"]]@jackstraw@overall.p.values
print(npcs <- max(which(p.values[,"Score"] <=0.05)))
npcs = 69
# Perform an integrated analysis
# Now we can run a single integrated analysis on all cells!
    
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(object) <- "integrated"

# Run the standard workflow for visualization and clustering
# Run the standard cca workflow for umap & tsne visualization
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = npcs, verbose = FALSE)
jpeg(paste0(path,"ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
print(ElbowPlot(object,ndims = 100))
dev.off()

object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))

object[["cca.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                             key = "ccaUMAP_", assay = DefaultAssay(object))
colnames(object[["cca.umap"]]@cell.embeddings) %<>% paste0("cca-",.)

object[["cca.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
                                             key = "cca-SNE_", assay = DefaultAssay(object))
colnames(object[["cca.tsne"]]@cell.embeddings) %<>% paste0("cca-",.)
object <- FindNeighbors(object,reduction = "umap",dims = 1:2, verbose = T)
object <- FindClusters(object, resolution = 0.8, algorithm= 1,  verbose = T, graph.name = "integrated_snn")

saveRDS(object, file = "data/mt_20220114.rds")


#======1.7 UMAP from raw pca =========================
format(object.size(object),unit = "GB")
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)

object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

#======1.8 UMAP from harmony =========================
DefaultAssay(object) = "SCT"

jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))

object[["harmony.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                 key = "harmonyUMAP_", assay = DefaultAssay(object))
colnames(object[["harmony.umap"]]@cell.embeddings) %<>% paste0("harmony-",.)

object[["harmony.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
                                                 key = "harmonytSNE_", assay = DefaultAssay(object))
colnames(object[["harmony.tsne"]]@cell.embeddings) %<>% paste0("harmony-",.)


object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))
object <- FindNeighbors(object,reduction = "umap",dims = 1:2, verbose = T)
object <- FindClusters(object, resolution = 0.8, algorithm= 1,  verbose = T)

UMAPPlot.1(object,group.by = "integrated_snn_res.0.8",do.print = T,label = T,label.repel = T)


saveRDS(object, file = "data/mt_20220114.rds")

