library(dplyr)
library(tidyr)
library(edgeR)
library(ggplot2)
library("RColorBrewer")

find_similar = function(gene_symb,pat_table,n = 20){
  key_gene = sym2ens(gene_symb)
  key_gene=gencode_clean(key_gene)
  if(is.na(key_gene)) stop("no such gene")
  key_pat = as.matrix(pat_table[key_gene,])
  key_pat = key_pat[1,]
  data = as.matrix(pat_table)
  cor_value=apply(data,1,function(x) cor(x,key_pat) )
  similar = pat_table[order(cor_value,decreasing=T),]
  genelist = rownames(similar[1:n,])
  return(genelist)
}

smad4_listbar = function(gene.expr){
  z = gather(gene.expr,libs,value,-gene)
  z = separate(z,libs,c("cellsdays","reps"),sep = -2) #D01A -> D01,A
  z = separate(z,cellsdays,c("cells","days"),sep = -4,remove = FALSE)
  
  p=ggplot(data = z,aes(x=cellsdays,y=value,fill= reps,shape = factor(cells)))+
    geom_bar(stat="identity",position="dodge") +
    facet_wrap( ~ gene,ncol = 5) #,scales = "free_y") #facet_wrap(~gene,scales="fixed")#
  q= p + theme(axis.text.x = element_text(angle = 45,hjust = 1))
  print(q)
}

smad4_single = function(gene.expr){
  z = gather(gene.expr,libs,value,-gene)
  z = separate(z,libs,c("cellsdays","reps"),sep = -2) #D01A -> D01,A
  z = separate(z,cellsdays,c("cells","days"),sep = -4,remove = FALSE)
  
  p=ggplot(data = z,aes(x=days,y=value,fill= reps,shape = factor(cells)))+
    geom_bar(stat="identity",position="dodge") +
    facet_wrap( ~ cells,ncol = 5) #,scales = "free_y") #facet_wrap(~gene,scales="fixed")#
  q= p + theme(axis.text.x = element_text(angle = 45,hjust = 1))
  print(q)
}



heatmapcolor=brewer.pal(n=3, name="RdYlBu")
#scale_fill_gradient2(low = heatmapcolor[3], mid = heatmapcolor[2],high = heatmapcolor[1],midpoint = 0.5)+

xgeom_heatmap  =   ggplot() +scale_fill_gradient(low = "brown1",high = "yellow") + 
  coord_fixed() + theme_classic()+labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())

geom_heatmap = ggplot()+
  scale_fill_gradient(low = "brown1",high = "yellow")+
  #scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0.5)+
  theme_classic()+labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())

lgeom_heatmap = ggplot()+
  scale_fill_gradient2(low = heatmapcolor[3], mid = heatmapcolor[2],high = heatmapcolor[1],midpoint = 0.5)+
  #scale_fill_gradient(low = "#FF6666",high = "#FFDB25") + 
  theme_minimal() + labs(x = NULL,y = "genes")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())

geom_heatmap_list = ggplot()+
  scale_fill_gradient2(low = heatmapcolor[3], mid = heatmapcolor[2],high = heatmapcolor[1],midpoint = 0.5)+
  coord_fixed() + theme_minimal()+labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())

smad4_heatmap = function(gene.expr) {
  mmat = gather(gene.expr,libs,value,-gene)
  p = ggplot()+
    geom_tile(data = mmat,aes(x = libs,y = gene, fill = value))+
    scale_fill_gradient(low = "brown1",high = "yellow")+
    #scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0.5)+
    scale_y_discrete(breaks = gene.expr$gene,label = gene.expr$gene,expand = c(0,0)) +
    coord_fixed() + theme_classic()+labs(x="",y="")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x = element_blank(),
          axis.ticks.y.left = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank())
  print(p)
}

smad4_heatmap_long = function(gene.expr) {
  z = gather(gene.expr,libs,value,-gene)
  z = separate(z,libs,c("cellsdays","reps"),sep = -2)
  z = separate(z,cellsdays,c("cells","days"),sep = -4,remove = FALSE)
  z = unite(z,daysreps,days,reps,sep = "",remove = F)
  
  p = lgeom_heatmap+
    geom_tile(data = z,aes(x = daysreps,y = gene, fill = value),color = "white") +
    facet_wrap(~cells)
  print(p)
}


smad4_heatmap_list = function(gene.expr) {
  z = gather(gene.expr,libs,value,-gene)
  z = separate(z,libs,c("cellsdays","reps"),sep = -2) #D01A -> D01,A
  z = separate(z,cellsdays,c("cells","days"),sep = -4,remove = FALSE)
  z = unite(z,daysreps,days,reps,sep = "",remove = F)
  
  p = geom_heatmap_list +
    scale_y_discrete(breaks = gene.expr$gene,label = gene.expr$gene,expand = c(0,0)) +
    geom_tile(data = z,aes(x = daysreps,y = gene, fill = value),color = "white")+
    facet_wrap(~cells)
  print(p)
}


gnorm = function(v,A,B){
  v.A = v[,A]
  v.B = v[,B]
  v.A = data.frame(t(apply(v.A,1,function(x) x/max(x))))
  v.B = data.frame(t(apply(v.B,1,function(x) x/max(x))))
  v[,A] = v.A
  v[,B] = v.B
  return(v)
}

load_smad4 = function(){
  file = "/Users/jiexu/kclab/SMAD4-RNA-seq/K.Chien_17_01-P8104/tables/smad4rnaseq_gene_tpm.txt"
  raw = read.delim(file,row.names = 1)
  keep = apply(raw,1,function(x) sum(x) > 0.5)
  rpkmtable = raw[keep,]
  return(rpkmtable)
}

load_smad4counts = function(){
  file = "/Users/jiexu/kclab/SMAD4-RNA-seq/K.Chien_17_01-P8104/tables/smad4_gene_merge.txt"
  #file = "/Users/jiexu/kclab/SMAD4-RNA-seq/K.Chien_17_01-P8104/tables/smad4rnaseq_gene_counts.txt"
  raw = read.delim(file,row.names = 1)
  keep = rowSums(cpm(raw)>1) >= 2
  counttable = raw[keep, ]
}

load_zic2cpm = function(){
  file = "/Users/jiexu/kclab/zic2-rna-seq/merged_gene_counts.txt"
  #file = "/Users/jiexu/kclab/SMAD4-RNA-seq/K.Chien_17_01-P8104/tables/smad4rnaseq_gene_counts.txt"
  raw = read.delim(file,row.names = 1)
  keep = rowSums(cpm(raw)>1) >= 2
  counttable = raw[keep, ]
  
  e = DGEList(counts=counttable)
  
  cpmtable = data.frame(cpm(e))
  
}

load_smad4clusters = function(){
  file = "/Users/jiexu/kclab/SMAD4-RNA-seq/K.Chien_17_01-P8104/my_results/cluster2/smad4_cluster_kmeans.txt"
  craw = read.delim(file,row.names = 26)
  keep = rownames(craw)
  raw = load_smad4counts()
  subcount = raw[keep,]
  A = rownames(craw)
  B = rownames(subcount)
  if(all(A==B))
    subcount$cluster = craw[,25]
  else
    stop("error in load_smad4cluster")
  return(subcount)
}

rearrange_smad4 = function(raw,A,B){
  result = raw
  libs = colnames(raw)
  A = grep('A$',libs)
  B = grep('B$',libs)
  L  = dim(raw)[2]
  result[,seq(1,L,2)] = raw[,A]
  result[,seq(2,L,2)] = raw[,B]
  rlibs = libs
  rlibs[seq(1,L,2)] = libs[A]
  rlibs[seq(2,L,2)] = libs[B]
  colnames(result) = rlibs
  return(result)
}


smad4_MDSedgeR=function(counttable){
  #this use edgeR
  e = DGEList(counts=counttable)
  mds = plotMDS(e)
  mds.points = data.frame(labels = rownames(mds$cmdscale.out),X1 = mds$x,X2 = mds$y)
  mds.points = separate(mds.points,labels,c("cellsdays","reps"),sep = -3)
  mds.points = separate(mds.points,cellsdays,c("cells","days"),sep = -4)
  
  cmd.plot = ggplot(data = mds.points,aes(x = X1,y = -X2)) +
    #geom_path(aes(group = cells),size = 2,alpha = I(1/4)) +
    geom_point(aes(shape = cells),color = "black",size = 7) +
    geom_point(aes(shape = cells, color = days),size = 5)
  print(cmd.plot)
}

smad4_mdsplot=function(show_data){
  cmat = cor(show_data)
  dmat = dist(t(cmat))#this gives the euclidean distance between different days
  #dmat = as.dist(1-cmat)
  cdmat = cmdscale(dmat,k = 4,eig = T,x.ret = T)
  #ic = isoMDS(dmat)
  #sa = sammon(dmat) 
  mds.points = data.frame(labels = rownames(cdmat$points),cdmat$points)
  mds.points = separate(mds.points,labels,c("cellsdays","reps"),sep = -3)
  mds.points = separate(mds.points,cellsdays,c("cells","days"),sep = -4)
  
  cmd.plot = ggplot(data = mds.points,aes(x = X1,y = -X2)) +
    #geom_path(aes(group = cells),size = 2,alpha = I(1/4)) +
    geom_point(aes(shape = cells),color = "black",size = 7) +
    geom_point(aes(shape = cells, color = days),size = 5)
  print(cmd.plot)
}

smad4_pcaplot=function(show_data){
  ir.pca = prcomp(show_data,center=T,scale.=F)
  pcs = data.frame(labels = rownames(ir.pca$rotation),ir.pca$rotation)
  pcs = separate(pcs,labels,c("cellsdays","reps"),sep = -3)
  pcs = separate(pcs,cellsdays,c("cells","days"),sep = -4)
  
  pcs.plot = ggplot(data = pcs,aes(x = PC1,y = PC2)) + 
    #geom_path(aes(group = cells),size = 2,alpha = I(1/4)) +
    geom_point(aes(shape = cells),color = "black",size = 7) +
    geom_point(aes(shape = cells, color = days),size = 5)
  print(pcs.plot)
}


smad4_cluster.pca=function(show_data){
  A = select(show_data,ends_with("A"))
  B = select(show_data,ends_with("B"))
  expr_data = (A+B)/2
  
  data=t(expr_data)
  ir.pca=prcomp(data,center=T,scale.=F)
  pcs=as.data.frame(ir.pca$rotation)
  pcs$cluster = show_data$cluster
  
  pc = ggplot()+ geom_point(data=pcs,aes(x=PC1,y=PC2,color=cluster),alpha=I(2/5))
  print(pc)
  print(pc+facet_wrap(~cluster,nrow = 2))
}


gg_smear = function(lrt,detags,title) {
  rns = rownames(lrt)
  yes = rns %in% detags
  p = ggplot() +
    geom_point(data = lrt[!yes,],aes(x = logCPM,y = logFC),alpha = I(1/5)) + 
    geom_point(data = lrt[detags,],aes(x = logCPM,y = logFC), #-PValue
               color = "red",alpha = I(1/3)) +
    labs(title = title)
  print(p)
}

gg_smear2 = function(lrt0,title,lfc = 3, pvalue = 1e-5) {
  lrt = lrt0$table
  rns = rownames(lrt)
  dt = decideTestsDGE(lrt0,lfc=lfc,p.value = pvalue)
  up_detags = rownames(lrt0$table)[dt > 0]
  down_detags = rownames(lrt0$table)[dt < 0]
  detags = union(up_detags,down_detags)
  yes = rns %in% detags
  p = ggplot(data = lrt,aes(x = logCPM,y = logFC,size = -log(PValue))) +
    geom_point(data = lrt[!yes,],alpha = I(1/7)) + 
    geom_point(data = lrt[up_detags,],color = "green",alpha = I(1/3)) +
    geom_point(data = lrt[down_detags,],color = "red",alpha = I(1/3)) +
    labs(title = title)+ theme_bw() #+ylim(-15,15)
  print(p)
}

gg_smear3 = function(lrt0,title,lfc = 3, pvalue = 1e-5) {
  lrt = lrt0$table
  rns = rownames(lrt)
  dt = decideTestsDGE(lrt0,lfc=lfc,p.value = pvalue)
  up_detags = rownames(lrt0$table)[dt > 0]
  down_detags = rownames(lrt0$table)[dt < 0]
  detags = union(up_detags,down_detags)
  yes = rns %in% detags
  p = ggplot(data = lrt,aes(x = logCPM,y = logFC,size = -log(PValue))) +
    geom_point(data = lrt[!yes,],alpha = I(1/7)) + 
    geom_point(data = lrt[up_detags,],color = "green",alpha = I(1/3)) +
    geom_text(data = lrt[up_detags,],aes(label = ens2sym(up_detags)),
              hjust = 0, nudge_x = 0.05,check_overlap = T,size = 3,color = 'green')+
    geom_point(data = lrt[down_detags,],color = "red",alpha = I(1/3)) +
    geom_text(data = lrt[down_detags,],aes(label = ens2sym(down_detags)),
              hjust = 0, nudge_x = 0.05,check_overlap = T,size = 3,color = 'red')+
    labs(title = title)+ theme_bw() + ylim(-15,15)
  print(p)
}

gg_volcano = function(ltable,dgenes){
  #dgenes = diff_genes(fit,my.contrasts[,'es03smad4'],'es03smad4',lfc = lfc,pvalue = pvalue)
  p = ggplot(data = ltable,aes(x = logFC,y = -log(PValue)))+
    geom_point()+
    geom_text(data = ltable[dgenes$up,],aes(label = ens2sym(gene)),
              hjust = 0, nudge_x = 0.05,check_overlap = T,size = 2,color = 'blue')+
    geom_text(data = ltable[dgenes$down,],aes(label = ens2sym(gene)),
              hjust = 0, nudge_x = 0.05,check_overlap = T,size = 2,color = 'red')
  print(p)
}

diff_genes = function(fit,contrast,title,lfc = 3, pvalue = 1e-5){
  lrt <- glmLRT(fit, contrast=contrast)
  dt = decideTestsDGE(lrt,lfc=lfc,p.value = pvalue)
  detags = rownames(lrt$table)[as.logical(dt)]
  #gg_smear(lrt,detags,title)
  #gg_smear2(lrt,title)
  gg_smear3(lrt,title)
  up = rownames(lrt$table)[dt == 1]
  down = rownames(lrt$table)[dt == -1]
  ans = list(up = up, down = down)
}



#this function is for reoder clusters for heatmap
smad4_sd=function(x){
  xm = 0.5*(x[seq(1,24,2)]+x[seq(2,24,2)])
  sd_0 = sd(xm[seq(1,12,4)])
  sd_1 = sd(xm[seq(2,12,4)])
  sd_2 = sd(xm[seq(3,12,4)])
  sd_3 = sd(xm[seq(4,12,4)])
  sds= c(sd_0,sd_1,sd_2,sd_3)
  return(max(sds))
}


#
write_list = function(filename,outlist){
  outf = file(description = paste0(filename),open = "wt")
  cat(paste(outlist,collapse="\n"),file=outf);
  cat("\n",file=outf)
  close(outf)
}

write_lists=function(filename,outlist){
  outf = file(description = paste0(filename),open = "wt")
  L = length(outlist)
  for(ii in seq(1,L)){
    templist = unlist(outlist[ii])
    cat(paste(templist,collapse="\t"),file=outf);
    cat("\n",file=outf)
  }
  close(outf)
}


