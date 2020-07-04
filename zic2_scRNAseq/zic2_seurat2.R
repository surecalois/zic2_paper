source('~/Dropbox/my_work_stuff/Resource/Scripts/some_r/my_rlibs/my_functions.R')
library(Seurat)
dropbox_dir = '/Users/jiexu/Dropbox'
#dropbox_dir ='/Users/xx/Dropbox'
annotation_file =paste0(dropbox_dir,'/my_work_stuff/Resource/expression_tables/zic2_scRNAseq/zic2_scRNAseq_annotation.txt')
molecules_file =paste0(dropbox_dir,'/my_work_stuff/Resource/expression_tables/zic2_scRNAseq/zic2_scRNAseq_clean.txt')

raw <- read.table(molecules_file, sep = "\t")
molecules = raw
x = ens2sym(row.names(molecules))
ix = !duplicated(ens2sym(x))
molecules = raw[ix,]
row.names(molecules) = ens2sym(row.names(molecules))

anno <- read.table(annotation_file, sep = "\t", header = TRUE)
row.names(anno) = anno$sample_id

sth = CreateSeuratObject(molecules,meta.data = anno)


mt_ribo_file = paste0(dropbox_dir,'/my_work_stuff/Resource/DNAinfo/list/Mito_ribosomal.txt')
s_ribo_file = paste0(dropbox_dir,'/my_work_stuff/Resource/DNAinfo/list/S_ribosomal.txt')
l_ribo_file = paste0(dropbox_dir,'/my_work_stuff/Resource/DNAinfo/list/L_ribosomal.txt')

mt_ribo = ens2sym(load_list(mt_ribo_file))
s_ribo = ens2sym(load_list(s_ribo_file))
l_ribo = ens2sym(load_list(l_ribo_file))


sth[["percent.mt"]] <- PercentageFeatureSet(sth, features = mt_ribo)

VlnPlot(sth, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)

plot1 <- FeatureScatter(sth, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sth, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


sth <- subset(sth, subset = nFeature_RNA > 6000 & nFeature_RNA < 12500)
sth <- NormalizeData(sth, normalization.method = "LogNormalize", scale.factor = 10000)
sth <- NormalizeData(sth)
sth <- FindVariableFeatures(sth, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(sth), 10)
plot1 <- VariableFeaturePlot(sth)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(sth)
sth <- ScaleData(sth, features = all.genes)

sth <- RunPCA(sth, features = VariableFeatures(object = sth))
PCAPlot(sth)
VizDimLoadings(sth, dims = 1:2, reduction = "pca")

DimPlot(sth, reduction = "pca")

sth <- RunTSNE(sth, features = VariableFeatures(object = sth))
TSNEPlot(sth)


DimHeatmap(sth, dims = 1, cells = 800, balanced = TRUE)

sth <- FindNeighbors(sth, dims = 1:10)
sth <- FindClusters(sth, resolution = 0.5)

sth <- RunUMAP(sth, dims = 1:10)
UMAPPlot(sth,shape.by = 'group')
UMAPPlot(sth,group.by = 'group',pt.size = 2)

DimPlot(sth, reduction = "umap")

FeaturePlot(sth, features = c('POU5F1','NANOG','SOX2','NODAL','MIXL1','T','EOMES','MESP1','TBX6'),
            shape.by = 'group',ncol = 3)

FeaturePlot(sth, features = c('MGMT','RP11-267L5.1','APLNR','LHX5','LIX1','RP3-467D16.3'),
            shape.by = 'group',ncol = 3)


sth.markers <- FindAllMarkers(sth, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sth.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- sth.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(sth, features = top10$gene)
