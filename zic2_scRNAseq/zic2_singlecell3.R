library(SingleCellExperiment)
library(scater)
library(scran)
library(SC3)
library(M3Drop)
options(stringsAsFactors = FALSE)

dropbox_dir = '/Users/jiexu/Dropbox'
#dropbox_dir ='/Users/xx/Dropbox'
annotation_file =paste0(dropbox_dir,'/my_work_stuff/Resource/expression_tables/zic2_scRNAseq/zic2_scRNAseq_annotation.txt')
molecules_file =paste0(dropbox_dir,'/my_work_stuff/Resource/expression_tables/zic2_scRNAseq/zic2_scRNAseq_clean.txt')
tsne_file = paste0(dropbox_dir,'/my_work_stuff/Resource/expression_tables/zic2_scRNAseq/zic2_tsne.txt')
#molecules_file = '/Users/jiexu/Desktop/working/zic2_2019/scRNA-seq/zic2_scRNAseq_clean.txt'
#annotation_file = '/Users/jiexu/Desktop/working/zic2_2019/scRNA-seq/zic2_annotation.txt'

#===============functions========================#
tsne_gene = function(zic2_mesp,gene = 'MESP1',tsne_name = 'tsne_fix'){
  zic2_tsne = data.frame(reducedDim(zic2_mesp,type = tsne_name))
  colnames(zic2_tsne) = c("tsne1","tsne2")
  zic2_tsne$group = zic2_mesp$group
  #zic2_tsne$gene = logcounts(zic2_mesp)[sym2ens(gene),]
  zic2_tsne$gene = logcounts(zic2_mesp)[gene,]
  midp = mean(zic2_tsne$gene[zic2_tsne$gene > 0])
  colnames(zic2_tsne)[4] = gene
  library(RColorBrewer)
  heatmapcolor=brewer.pal(n=3, name="RdYlBu")
  p = ggplot(zic2_tsne,aes_string(x = "tsne1",y = "tsne2",shape = "group",color = gene))+geom_point()+
    scale_color_gradient2(low = heatmapcolor[3], mid = heatmapcolor[2],high = heatmapcolor[1],
                          midpoint = midp) + 
    theme_classic()
  return(p)
}

tsne_genes = function(zic2_mesp,genes,tsne_name = 'tsne_fix'){
  zic2_tsne = data.frame(reducedDim(zic2_mesp,type = tsne_name))
  colnames(zic2_tsne) = c("tsne1","tsne2")
  zic2_tsne$group = zic2_mesp$group
  #inlist = sym2ens(genes)
  inlist = genes
  keep.name = match(inlist,rownames(zic2_mesp)) 
  genes = genes[!is.na(keep.name)]
  genes = sub('-','.',genes)
  keep.name = keep.name[!is.na(keep.name)]
  w = logcounts(zic2_mesp)[keep.name,]#you need to change sym2ens to other 
  #midp = mean(zic2_tsne$gene[zic2_tsne$gene > 0])
  w = norm.by.max(w) #map max cell to value 1
  w = t(w)
  colnames(w) = genes
  zic2_tsne = data.frame(zic2_tsne,w)
  zz = gather(zic2_tsne,key = "gene",value = "value",3+seq(1,length(genes)))
  zz$gene = factor(zz$gene,levels = genes)
  library(RColorBrewer)
  heatmapcolor=brewer.pal(n=3, name="RdYlBu")
  p = ggplot(zz,aes_string(x = "tsne1",y = "tsne2",shape = "group",color = "value"))+geom_point()+
    scale_color_gradient2(low = heatmapcolor[3], mid = heatmapcolor[2],high = heatmapcolor[1],
                          midpoint = 0.5) + 
    theme_classic()+
    theme(strip.background = element_rect(colour="#FFFFFF", fill="#FFFFFF"))+ 
    facet_wrap(~gene) #,ncol = 1
  
  return(p)
}

#===============end of functions=================#


molecules <- read.table(molecules_file, sep = "\t")
anno <- read.table(annotation_file, sep = "\t", header = TRUE)

zic2_mesp = SingleCellExperiment(
  assays = list(counts = as.matrix(molecules)), 
  colData = anno
)

keep_feature <- rowSums(counts(zic2_mesp) > 0) > 0
zic2_mesp <- zic2_mesp[keep_feature, ]

mt_ribo_file = paste0(dropbox_dir,'/my_work_stuff/Resource/DNAinfo/list/Mito_ribosomal.txt')
s_ribo_file = paste0(dropbox_dir,'/my_work_stuff/Resource/DNAinfo/list/S_ribosomal.txt')
l_ribo_file = paste0(dropbox_dir,'/my_work_stuff/Resource/DNAinfo/list/L_ribosomal.txt')

mt_ribo = sym2ens(load_list(mt_ribo_file)) 
s_ribo = sym2ens(load_list(s_ribo_file)) 
l_ribo = sym2ens(load_list(l_ribo_file)) 

isSpike(zic2_mesp, "mt_ribosomal") <- rownames(zic2_mesp) %in% mt_ribo
isSpike(zic2_mesp, "s_ribosomal") <- rownames(zic2_mesp) %in% s_ribo
isSpike(zic2_mesp, "l_ribosomal") <- rownames(zic2_mesp) %in% l_ribo
#check it in the zic2_mesp in the spikeNames tag


zic2_mesp = calculateQCMetrics(zic2_mesp,use_spikes = T)
#names(colData(zic2_mesp)) # this is for the cells
#names(rowData(zic2_mesp)) # this is for the genes


#the distribution of the total counts of each cell
hist(
  zic2_mesp$total_counts,
  breaks = 100 # why choose 100
)
abline(v = .12e6, col = "red")
filter_by_total_counts <- (zic2_mesp$total_counts > .12e6 & zic2_mesp$total_counts < 1e6)

#the distribution of the number of (features = genes) detected in each cell
hist(
  zic2_mesp$total_features_by_counts,
  breaks = 100 # why choose 100
)
abline(v = 7e3, col = "red")
filter_by_expr_features <- (zic2_mesp$total_features_by_counts > 7e3)

zic2_mesp$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts
)
#table(zic2_mesp$group[zic2_mesp$use]) #to show many cells in the dataset


filter_genes <- apply(
  counts(zic2_mesp), 
  1, 
  function(x) length(x[x > 1]) >= 5
)
rowData(zic2_mesp)$use <- filter_genes
#table(filter_genes) # to display how many pass and how many fail


zic2_mesp.qc <- zic2_mesp[rowData(zic2_mesp)$use, colData(zic2_mesp)$use]
zic2_mesp.bk = zic2_mesp
zic2_mesp = zic2_mesp.qc


#===========normalized with scran=======================#
clusters = quickCluster(zic2_mesp, min.mean=0.1, method="igraph")
zic2_mesp      = computeSumFactors(zic2_mesp, cluster=clusters, min.mean=0.1)
zic2_mesp <- computeSpikeFactors(zic2_mesp, general.use=FALSE)

summary(sizeFactors(zic2_mesp))
sizeFactorNames(zic2_mesp)
summary(sizeFactors(zic2_mesp,"mt_ribosomal"))
summary(sizeFactors(zic2_mesp,"s_ribosomal"))
summary(sizeFactors(zic2_mesp,"l_ribosomal"))

plot(zic2_mesp$total_counts, sizeFactors(zic2_mesp), log="xy", 
     xlab="total counts", ylab="size factors", 
     cex=0.3, pch=20, col=rgb(0.1,0.2,0.7,0.3))

assayNames(zic2_mesp)
zic2_mesp = normalize(zic2_mesp) #log(count/size.factor + 1)
zic2_mesp = normalize(zic2_mesp,return_log = F) #then this calculate "normcounts"
assay(zic2_mesp, "cpm")<- calculateCPM(zic2_mesp)
assay(zic2_mesp, "logcounts_raw") <- log2(counts(zic2_mesp) + 1)
assay(zic2_mesp, "logcounts_cal") <- log2(normcounts(zic2_mesp) + 1) #this is the as normalize
assay(zic2_mesp, "logcounts_cpm") <- log2(calculateCPM(zic2_mesp) + 1)
assayNames(zic2_mesp)

#===========load the tsne display you stored before from a file=============#
#zic2_mesp = runTSNE(zic2_mesp)
#reducedDimNames(zic2_mesp)
tsne_fix = read.table(tsne_file, sep = "\t",header = T)
tsne_m = as.matrix(tsne_fix[,2:3])
rownames(tsne_m) = tsne_fix[,1]
reducedDim(zic2_mesp,"tsne_fix") = tsne_m
reducedDimNames(zic2_mesp)

#==============clustering==================#
rownames(zic2_mesp) = ens2sym(rownames(zic2_mesp))
rowData(zic2_mesp)$feature_symbol = rownames(zic2_mesp)

zic2_mesp <- sc3_estimate_k(zic2_mesp)
metadata(zic2_mesp)$sc3$k_estimation

sc3_k_range = 8:10
zic2_mesp <- sc3(zic2_mesp, ks = sc3_k_range, biology = TRUE) #this take very long time (3~5 mins)

sc3_k = 9
#for(sc3_k in sc3_k_range){
 # sc3_plot_consensus(zic2_mesp, k = sc3_k, show_pdata = "group")
  #sc3_plot_expression(zic2_mesp, k = sc3_k, show_pdata = "group")
  sc3_plot_markers(zic2_mesp, k = sc3_k, show_pdata = c("group"))
  #sc3_plot_de_genes(zic2_mesp,k = the_k, show_pdata = "group")
#}
  
  sc3_clusters_str = paste0("sc3_",sc3_k,"_clusters")
  p = plotReducedDim(zic2_mesp, "tsne_fix", colour_by = sc3_clusters_str,size_by = "total_counts")#,shape_by = "group")
  p =p +  scale_fill_manual(values=c("#ccaea7","#ccaea7","#ccaea7", "#02d9df", "#f6aa00","#b0abff","#ff8aff","#00db98","#00cbff"))
  print(p)
  
  p = plotReducedDim(zic2_mesp, "tsne_fix", colour_by = "group",size_by = "total_counts")
  p = p + scale_fill_manual(values=c("#bec100", "#ff9189", "#5cd131"))
  print(p)
  
tsne_gene(zic2_mesp,'EOMES')


tsne_genes(zic2_mesp,c('POU5F1','NANOG','SOX2','NODAL','MIXL1','T','EOMES','MESP1','TBX6'))


x = rowData(zic2_mesp)$sc3_10_markers_padj
ox = order(x)
y = ens2sym(rownames(zic2_mesp)[ox[1:10]])

tsne_genes(zic2_mesp,y)

markerts = c('MGMT','RP11-267L5.1','APLNR','LHX5','LIX1','RP3-467D16.3')
  tsne_genes(zic2_mesp,markerts)

row_data <- rowData(zic2_mesp)
zic2_sc3 = row_data[ , grep("sc3_", colnames(row_data))]
zic2_sc3 = zic2_sc3[zic2_sc3$sc3_gene_filter,]
zic2_sc3 = zic2_sc3[,-1]


names(colData(zic2_mesp)) #about the cells
names(rowData(zic2_mesp)) #about the genes
