library(dplyr)
library(tidyr)
library(edgeR)
library(ggplot2)
library(RColorBrewer)
heatmapcolor=brewer.pal(n=3, name="RdYlBu")

usr_home = path.expand('~')

gencode_dir = paste0(usr_home,"/Dropbox/my_work_stuff/Resource/DNAinfo/gencode")
gene_map_v24 = paste0(gencode_dir,"/V24/gene_map_v24.txt")
ens2gencode_v24 = paste0(gencode_dir,"/V24/ens2gencode_v24.txt")

gene_map_v26 = paste0(gencode_dir,"/V26/gene_map_v26lift37.txt")
ens2gencode_v26 = paste0(gencode_dir,"/V26/ens2gencode_v26lift37.txt")

#============================begin of filters============================================#
#keep = rowSums(cpm(raw)>1) >= 2
#keep = apply(raw,1,function(x) sum(x) > 0.5)
# why don't we use some genes as reference

filter.low = function(expr.table,th = 1){
  keep = apply(expr.table,1,function(x) max2(x) > th)
}

filter.dyn = function(expr.table, th = 0.1){ #there isn't really a good way to do this
  keep = apply(expr.table,1,function(x) sd(x)/(1+max(x)) > th)
}

filter.sd = function(expr.table, th = 0.75){ #there isn't really a good way to do this
  sd.vector = apply(expr.table,1,sd)
  th.value = quantile(sd.vector,probs = th)
  keep = sd.vector > th.value
}

filter.similar = function(expr.table,th = 0.8){
  keep = apply(expr.table,1,function(x) ps2(x) > th)
}

# tpm.max = apply(a,1,max)
# tpm.sd = apply(a,1,sd)
# 
# tpm.key = tpm.sd/(1+tpm.max) # why you choose 1, not 0.01
# 
# ht = hist(tpm.key,breaks = 1000)
# plot(density(tpm.key))

##some row function for data.frame using apply
#these function require A,B,A,B sequence
max2 = function(x){
  n = length(x)
  max.a = max(x[seq(1,n,2)])
  max.b = max(x[seq(2,n,2)])
  return(min(max.a,max.b))
}

ps2 = function(x) {
  n = length(x)
  if(sd2(x) > 2)
    ans = cor(x[seq(1,n,2)],x[seq(2,n,2)])
  else
    ans = 0
  return(ans) #measure the repeatability of the pattern. 
}

#measure the dynamic of the pattern. 
sd2 = function(x) {
  n = length(x)
  return(sqrt(sd(x[seq(1,n,2)])*sd(x[seq(2,n,2)])))
}

#============================end of filters============================================#

#============================begin of cluster============================================#

kmeans_it = function(expr.table,cluster_num = 12,plot.kmc = T){
  # I think you should do the filtering before feed in this function
  # a = expr.table
  # a = a[filter.low(a,th1) ,]
  # a = a[filter.similar(a, th2) ,]
  # a = a[filter.dyn(a, th3) ,]
  # expr.table = a
  
  pat_table = norm.by.group(expr.table)
  vv = pat_table[complete.cases(pat_table),] #remove NAs rows
  
  #set.seed(12345)
  km = kmeans(vv,cluster_num, iter.max = 500000)
  
  kmc = t(apply(km$centers,1,function(x) x/sum(x)))
  wmi = apply(kmc,1,sd2)
  wi = order(wmi,decreasing=F) #on heatmap, the first one goes to bottom
  
  if(plot.kmc) {
    kmc.df = data.frame(km$centers)
    kmc.df$gene = factor(rownames(kmc.df),levels = wi, ordered = T)
    heatmap_list(kmc.df)
  }
  
  
  svv = data.frame(vv,cluster = factor(km$cluster,levels = wi,ordered = T))
  
 # print(svv[1:5,20:22])
 # plot(as.numeric(svv[1,1:20]))
#  plot(kmc[svv[1,22],])
  dct = function(x) {
    y = rev(x)
    ct = y[1]
   # if(ct %in% wi){
      xt = as.numeric(y[-1])
      xt = rev(xt)
      yt = as.numeric(kmc[ct,])
      ans = cor(xt,yt)
  #  }else{
  #    ans = 0
  #  }
    return(ans)
  }
  svv$sp = apply(svv,1,dct)
  svv = svv[order(svv$cluster),]
  #svv = svv[order(svv$cluster,-svv$sp),]
  #svv = svv[svv$sp > .70,]
}

show.2d = function(mds.points){
  cmd.plot = ggplot(data = mds.points,aes(x = mds.points[,5],y = mds.points[,6])) +
    geom_point(aes(shape = cells, color = days,size = 5)) #+ 
   # scale_shape_manual(values = levels(factor(mds.points$cells)))
  print(cmd.plot)
}

MDS.others = function(counttable){
  cmat = cor(counttable)
  dmat = dist(t(cmat))#this gives the euclidean distance between different days
  #dmat = as.dist(1-cmat)
  cdmat = cmdscale(dmat,k = 2,eig = T,x.ret = T)
  #ic = isoMDS(dmat)
  #sa = sammon(dmat) 
  mds.points = data.frame(libs = rownames(cdmat$points),cdmat$points)
  mds.points = separate(mds.points,libs,c("cells","days","reps"),sep ="\\.",remove = FALSE)
}

MDS.edgeR=function(counttable){
  #this use edgeR
  e = DGEList(counts=counttable)
  mds = plotMDS(e)
  mds.points = data.frame(libs = rownames(mds$cmdscale.out),X1 = mds$x,X2 = mds$y)
  mds.points = separate(mds.points,libs,c("cells","days","reps"),sep ="\\.",remove = FALSE)
}

PCA.columns=function(show_data){
  ir.pca = prcomp(show_data,center=T,scale.=F)
  pcs = data.frame(libs = rownames(ir.pca$rotation),ir.pca$rotation)
  pcs = separate(pcs,libs,c("cells","days","reps"),sep ="\\.",remove = FALSE)
}

PCA.gene=function(expr.table){
  A = select(expr.table,ends_with("A"))
  B = select(expr.table,ends_with("B"))
  show_data = (A+B)/2 #there is so problem abount this average.
  
  day_label = unlist(strsplit(colnames(show_data),"\\."))
  max.index = apply(show_data,1,which.max)
  max.value = apply(show_data,1,max)
  
  data=t(show_data)
  ir.pca=prcomp(data,center=T,scale.=T)
  pcs=as.data.frame(ir.pca$rotation)
  pcs$day=day_label[max.index]
  pcs$maxv=log10(max.value)
  
  pc = ggplot()+ geom_point(data=pcs,aes(x=PC1,y=PC2,size=maxv,color=day),alpha=I(2/5))
  print(pc)
  print(pc+facet_wrap(~day,nrow = 2))
  
  pcasdev = data.frame(sdev = ir.pca$sdev)
  sdns = paste0("PC",seq(1,10))
  pcasdev$component = factor(sdns,levels = sdns,ordered = T)
  p = ggplot()+
    geom_bar(data = pcasdev,aes(x = component,y = sdev),
             fill = "lightblue",stat = "identity",width = 0.7)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
}


#============================end of cluster==============================================#




#============================begin of DEA================================================#
make.design = function(table.colnames){
  gp = table.colnames
  w = unlist(strsplit(gp,"\\."))
  wl = length(w)
  
  my.cells = c('es03','smad4','zic2','isl1','nkx25','niwp')
  cell.regex = paste(my.cells,collapse = '|')
  cell = w[grep(cell.regex,w)]
  #day = w[grep('^d\\d+$',w)]
  day = w[grep('^\\D\\d+$',w)]
  rep = w[grep('[AB]',w)]
  
  groups = paste(cell,day,sep = ".")
  libs = gp
  sapply(list(libs,cell,day,rep,groups),length)
  expriment = data.frame(libs = libs, cells = cell, days = day, reps = rep, groups = groups)
  
  design = model.matrix(~0+groups,data = expriment)
  groups = factor(groups)
  colnames(design)=levels(groups)
  return(design)
}

make.groups = function(table.colnames){
  
  gp = table.colnames
  w = unlist(strsplit(gp,"\\."))
  wl = length(w)
    
  my.cells = c('es03','smad4','zic2','isl1','nkx25','niwp')
  cell.regex = paste(my.cells,collapse = '|')
  cell = w[grep(cell.regex,w)]
  day = w[grep('^d\\d+$',w)]
  rep = w[grep('[AB]',w)]
  groups = paste(cell,day,sep = ".")
  
}

make.gmlFit = function(counttable){
  gp = names(counttable)
  design = make.design(gp)
  groups = make.groups(gp)
  
  e = DGEList(counts=counttable,group = groups)
  e = calcNormFactors(e)
  e = estimateGLMCommonDisp(e, design)
  e = estimateGLMTrendedDisp(e,design)
  e = estimateGLMTagwiseDisp(e,design)
  fit = glmFit(e, design)
}


diff_genes = function(fit.result,contrast,some.genes,lfc = 3, pvalue = 1e-5,plot.md = T){
#example of contrast: 
#  my.contrast = makeContrasts(somename = es03.d0 - (zic2c1.d0+zic2c2.d0)/2, levels = design)
#  use columns of the contrasts: my.contrast[,1] or my.contrast[,somename]
  
  lrt <- glmLRT(fit.result, contrast=contrast)
  if(plot.md) diff_mdplot(lrt,some.genes,title = colnames(z.contrast))
  dt = decideTestsDGE(lrt,lfc=lfc,p.value = pvalue)
  detags = rownames(lrt$table)[as.logical(dt)]
  up = rownames(lrt$table)[dt == 1]
  down = rownames(lrt$table)[dt == -1]
  ans = list(up = up, down = down)
  return(ans)
}

diff_table = function(fit.result,contrast,lfc = 3, pvalue = 1e-5){
  #example of contrast: 
  #  my.contrast = makeContrasts(somename = es03.d0 - (zic2c1.d0+zic2c2.d0)/2, levels = design)
  #  use columns of the contrasts: my.contrast[,1] or my.contrast[,somename]
  
  lrt <- glmLRT(fit.result, contrast=contrast)
  table = lrt$table
  table$ens = rownames(table)
  table$FDR = p.adjust(table$PValue)
  result = decideTestsDGE(lrt,lfc=lfc,p.value = pvalue)
  table$result = as.vector(result)
  return(table)
}

diff_tops = function(diff.table,up = 15, down = 15){
  up.table = diff.table[diff.table$result == 1,]
  up.top = top_n(up.table,-up,PValue) # -up selection the bottom, smallest PValue
  
  down.table = diff.table[diff.table$result == -1,]
  down.top = top_n(down.table,-down,PValue)
  tops = rbind(up.top,down.top)
  rownames(tops) = tops$ens
  tops$gene = ens2sym(tops$ens)
  tops$ens = c()
  return(tops)
}


diff_mdplot = function(lrt0,some.genes,title = "MD.plot",lfc = 3, pvalue = 1e-5) {
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
    labs(title = title)+ theme_bw() + ylim(-15,15)
  
  if(!missing(some.genes)){
    p = p + geom_text(data = lrt[some.genes,],aes(label = ens2sym(gencode_clean(some.genes))),
              hjust = 0, nudge_x = 0.05,check_overlap = T,size = 3,color = 'blue')
  }
  
  print(p)
}

#==============================end of DES================================================#

#============================begin of list bar_plot=====================================#

giwi_listbar = function(gene.expr,nc = 5){
  zz = gather(gene.expr,libs,value,-gene)
  z = separate(zz,libs,c("cells","days","reps")) #D01A -> D01,A
  
  p=ggplot(data = z,aes(x=days,y=value,fill=reps))+geom_bar(stat="identity",position="dodge") +
    facet_wrap( ~ gene,ncol = nc)#,scales = "free_y" )#facet_wrap(~gene,scales="fixed")#
  q= p + theme(axis.text.x = element_text(angle = 45,hjust = 1))
  print(q)
}

#============================end of list bar_plot=====================================#


#============================begin of heatmaps===========================================#
#define some ggplot heatmap thing
#usage gg_heatmap+geom_tile(data = mmat,aes(x = libs,y = gene, fill = value))

geom_heatmap_long = ggplot()+
  scale_fill_gradient2(low = heatmapcolor[3], mid = heatmapcolor[2],high = heatmapcolor[1],midpoint = 0.5)+
  theme_minimal() + labs(x = NULL,y = "genes")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid = element_blank())

geom_heatmap_list = ggplot()+
  scale_fill_gradient2(low = heatmapcolor[3], mid = heatmapcolor[2],high = heatmapcolor[1],midpoint = 0.5)+
  #coord_fixed() + 
  theme_minimal()+labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid = element_blank())

#now define some easy-to-use function
prep_table = function(expr.table,check.complete = F){
  if(check.complete){
    z = complete.cases(expr.table)
    if(!all(z)) cat("some genes were not in the table --JX\n")
    expr.table = expr.table[z,]
  }
  rnames = rownames(expr.table)
  genes = ens2sym(rnames)
  dc = duplicated(genes)
  dropnames = rnames[dc]
  if(length(dropnames) > 10 ){
    k = 10
  }else{
    k = length(dropnames)
  }
  if(k > 0) cat("I drop these genes from the table,weird stuff: --JX\n",paste(dropnames[1:k],collapse = '\n'),"...")
  genes = genes[!dc]
  expr.table = expr.table[!dc,]
  expr.table$gene =  factor(genes,level = genes, ordered = T)
  return(expr.table)
}

heatmap_list = function(gene.expr,print.p = T) {
  if(!('gene' %in% names(gene.expr))){
    genes = rownames(gene.expr)
    gene.expr$gene = factor(genes,level = genes, ordered = T)
    warning("you didn't prep the table, so I will use rownames.")
  }
  mmat = gather(gene.expr,libs,value,-gene)
  p = geom_heatmap_list + 
    geom_tile(data = mmat,aes(x = libs,y = gene, fill = value),color = 'white') + 
    scale_y_discrete(breaks = gene.expr$gene,label = gene.expr$gene,expand = c(0,0))
  if(print.p) {
    print(p)
  }else{
    return(p)
  }
}

heatmap_long = function(gene.expr,print.p = T) {
  if(!('gene' %in% names(gene.expr))){
    genes = rownames(gene.expr)
    gene.expr$gene = factor(genes,level = genes, ordered = T)
    warning("you didn't prep the table, so I will use rownames.")
  }
  mmat = gather(gene.expr,libs,value,-gene)
  p = geom_heatmap_long + 
    geom_tile(data = mmat,aes(x = libs,y = gene, fill = value)) + 
    scale_y_discrete(breaks = gene.expr$gene,label = gene.expr$gene,expand = c(0,0))
  if(print.p) {
    print(p)
  }else{
    return(p)
  }
}

wrap_heatmap_list = function(gene.expr,groups,print.p = T) {
  if(!('gene' %in% names(gene.expr))){
    genes = rownames(gene.expr)
    gene.expr$gene = factor(genes,level = genes, ordered = T)
    warning("you didn't prep the table, so I will use rownames.")
  }
  z = gather(gene.expr,libs,value,-gene)
  z$groups = 'others'
  for(g in groups){
    z$groups[grep(g,z$libs)] = g
  }
  z = separate(z,libs,c("cells","days","reps"),sep ="\\.",remove = FALSE)
  #z = unite(z,daysreps,days,reps,sep = ".",remove = F)
  
  p = geom_heatmap_list +
    scale_y_discrete(breaks = gene.expr$gene,label = gene.expr$gene,expand = c(0,0)) +
    geom_tile(data = z,aes(x = libs,y = gene, fill = value),color = "white")+
    facet_wrap(~groups,nrow = 1, scales = "free_x")
  if(print.p) {
    print(p)
  }else{
    return(p)
  }
}

wrap_heatmap_long = function(gene.expr,groups,print.p = T) {
  if(!('gene' %in% names(gene.expr))){
    genes = rownames(gene.expr)
    gene.expr$gene = factor(genes,level = genes, ordered = T)
    warning("you didn't prep the table, so I will use rownames. --JX")
  }
  z = gather(gene.expr,libs,value,-gene)
  z$groups = 'others'
  for(g in groups){
    z$groups[grep(g,z$libs)] = g
  }
  z = separate(z,libs,c("cells","days","reps"),sep ="\\.",remove = FALSE)
  #z = unite(z,daysreps,days,reps,sep = ".",remove = F)
  
  p = geom_heatmap_long +
    scale_y_discrete(breaks = gene.expr$gene,label = gene.expr$gene,expand = c(0,0)) +
    geom_tile(data = z,aes(x = libs,y = gene, fill = value))+
    facet_wrap(~groups,nrow = 1, scales = "free_x")
  if(print.p) {
    print(p)
  }else{
    return(p)
  }
}


#============================end of heatmaps=============================================#

#============================begin of table normalized===================================#
norm.by.sum = function(some.df){
  out=data.frame(t(apply(some.df,1,function(x) x/sum(x))))
}

norm.by.max = function(some.df){
  out=data.frame(t(apply(some.df,1,function(x) x/max(x))))
}

norm.by.mean = function(some.df){
  out=data.frame(t(apply(some.df,1,function(x) x/mean(x))))
}

norm.by.median = function(some.df){
  out=data.frame(t(apply(some.df,1,function(x) x/median(x))))
}

# A = grep('A$',names(b))
# B = grep('B$',names(b))
norm.by.group = function(v,A,B,f = norm.by.max){
  if(missing(A)) A = grep('A$',names(v))
  if(missing(B)) B = grep('B$',names(v))
  v.A = v[,A]
  v.B = v[,B]
  v.A = f(v.A)
  v.B = f(v.B)
  v[,A] = v.A
  v[,B] = v.B
  return(v)
}

#============================end of table normalized===================================#



#=============================begin of list generators===================================#

load_list = function(file){
  
  if(missing(file)) file = paste0(usr_home,"/Desktop/templist.txt")
  
  #file = "/Users/jiexu/genome/annotation/list/FGFs.txt"
  #file = "/Users/jiexu/genome/annotation/list/tgfbeta.txt"
  #file = "/Users/jiexu/genome/annotation/list/wnt.txt"
  if(file == 'CMD') file = paste0(usr_home,
                                  "/Dropbox/my_work_stuff/Resource/DNAinfo/list/my_cardiacDiff.txt")
  
  listfile <- read.table(file, header=F,as.is = T);
  inlist = listfile$V1
  genelist=sym2ens(inlist); #you can also load ensembl ids now
}

find_similar = function(gene_symb,pat_table,n = 20){
  key_gene = sym2ens(gene_symb)
  #key_gene=gencode_clean(key_gene)
  if(is.na(key_gene)) stop("no such gene")
  key_pat = as.matrix(pat_table[key_gene,])
  key_pat = key_pat[1,]
  data = as.matrix(pat_table)
  cor_value=apply(data,1,function(x) cor(x,key_pat) )
  similar = pat_table[order(cor_value,decreasing=T),]
  genelist = rownames(similar[1:n,])
  return(genelist)
}


#============================end of list generators=============================#



#=============================begin of symbol converting system=================#


if(!exists("map24")){
  map24 = read.table(gene_map_v24)
  names(map24) = c("ens","sym","type")
  rownames(map24) = map24$ens
}

if(!exists("map26")){
  map26 = read.table(gene_map_v26)
  names(map26) = c("ens","sym","type")
  rownames(map26) = map26$ens
}

if(!exists("map24ens")){
  map24ens = read.table(ens2gencode_v24)
  names(map24ens) = c("ens","gencode","sym")
  rownames(map24ens) = map24ens$ens
}

if(!exists("map26ens")){
  map26ens = read.table(ens2gencode_v26)
  names(map26ens) = c("ens","gencode","sym","type")
  rownames(map26ens) = map26ens$ens
}

map24 = map26 #a stupid way to upgrade
map24ens = map26ens #



gencode2sym = function(ens_list,mapfile = map24){
  sym_list = mapfile[as.character(ens_list),"sym"]
  sym_list = as.character(sym_list)
  ind = is.na(sym_list)
  sym_list[ind] = as.character(ens_list)[ind]
  return(sym_list)
}

sym2gencode = function(sym_list,mapfile = map24){
  # hit = mapfile$sym %in% as.character(sym_list)
  hit = match(as.character(sym_list),map24$sym)
  ens_list = mapfile[hit,"ens"]
  return(as.character(ens_list))
}


ens2gencode = function(ens_list,mapfile = map24ens){
  sym_list = mapfile[as.character(ens_list),"gencode"]
  return(as.character(sym_list))
}

ens2sym = function(ens_list,mapfile = map24ens){
  if(nchar(ens_list[1]) > 15) ens_list = substr(ens_list,1,15) #now it can be used in both cases.
  sym_list = mapfile[as.character(ens_list),"sym"]
  sym_list = as.character(sym_list)
  ind = is.na(sym_list)
  sym_list[ind] = as.character(ens_list)[ind]
  return(sym_list)
}

sym2ens = function(sym_list,mapfile = map24){
  # hit = mapfile$sym %in% as.character(sym_list)
  ens_in = grep("^ENSG\\d{11}",sym_list)
  hit = match(as.character(sym_list),map24$sym)
  ens_list = mapfile[hit,"ens"]
  temp = as.character(ens_list)
  temp[ens_in] = sym_list[ens_in]
  return(substr(temp,1,15))
}

gencode_clean = function(ens_list){
  the_list = c()
  for(ens in ens_list){
    the_list = c(the_list,substr(ens,1,15))
  }
  return(the_list)
}

#======================end of symbol converting system===================#


#============================begin of data file location===================#
gencode_dir = paste0(usr_home,"/Dropbox/my_work_stuff/Resource/DNAinfo/gencode")

dir = paste0(usr_home,"/Dropbox/my_work_stuff/Resource/expression_tables/core/")
screen_file = paste0(dir,"core-zlib2all.txt")

giwi_tpm_file = paste0(dir,"giwi2_seq_st_tpm.txt")
giwi_fpkm_file = paste0(dir,"giwi2_seq_st_fpkm.txt")
giwi_count_file =paste0(dir,"giwi2_seq_fc_counts.txt")

mut_tpm_file = paste0(dir,"mutant_seq_st_tpm.txt")
mut_fpkm_file = paste0(dir,"mutant_seq_st_fpkm.txt")
mut_count_file =paste0(dir,"mutant_seq_fc_counts.txt")

msc_count_file = paste0(dir,"day6_msc.txt")

#=============================end of data file location======================#



#======================begin of data loading functions===================#

load_screen = function(){
  sgs=c("isl1","mesp1")
  gps=c("neg","pos")
  rpn=c("A","B","C")
  gp_rp = paste(gl(2,1,length = 6,label = gps),
                gl(3,2,length = 6,label = rpn),sep = '.')
  isl1.lib=paste("isl1",gp_rp,sep=".")
  mesp1.lib=paste("mesp1",gp_rp,sep=".")
  hESC.lib=paste("zlib2",rpn,sep=".")
  
  #datafile = "/Users/jiexu/Dropbox/PhD/Resource/tabels/screen/core-zlib2all.txt"
  datafile = screen_file
  raw <- read.table(datafile, header=F);
  names(raw)=c("chr","st","ed","sgRNA","gene","sig",isl1.lib,mesp1.lib,hESC.lib,"total","pam","seq")
  rownames(raw)=raw[,"sgRNA"]
  return(raw)
}

load_giwi2 = function(){
  print("giwi2 loaded")
  datafile = giwi_tpm_file
  #datafile = giwi_fpkm_file
  #datafile="/Users/jiexu/Desktop/owncloud/kclab/expression_tables/core/giwi2_cf_gene.txt"
  raw = read.table(datafile, header=T,row.names = 1);
  old_names = names(raw)
  days = substr(old_names,1,3)
  reps = substr(old_names,4,4)
  names(raw) = paste('es03',days,reps,sep = '.')
  keep = apply(raw,1,function(x) sum(x) > 0.5)
  rpkmtable = raw[keep,]
  rownames(rpkmtable) = substr(rownames(rpkmtable),1,15)
  return(rpkmtable)
}

load_giwi2counts = function(){
  datafile = giwi_count_file
  #datafile="/Users/jiexu/Desktop/owncloud/kclab/expression_tables/core/giwi2_fc_counts.txt"
  raw = read.table(datafile, header=T,row.names = 1)
  old_names = names(raw)
  days = substr(old_names,1,3)
  reps = substr(old_names,4,4)
  names(raw) = paste('es03',days,reps,sep = '.')
  keep = rowSums(edgeR::cpm	(raw)>1) >= 2
  counttable = raw[keep, ]
  rownames(counttable) = substr(rownames(counttable),1,15)
  return(counttable)
}

load_giwi2cpm = function(){
  file = giwi_count_file
  raw = read.delim(file,row.names = 1)
  old_names = names(raw)
  days = substr(old_names,1,3)
  reps = substr(old_names,4,4)
  names(raw) = paste('es03',days,reps,sep = '.')
  keep = rowSums(edgeR::cpm	(raw)>1) >= 2
  counttable = raw[keep, ]
  e = DGEList(counts=counttable)
  cpmtable = data.frame(edgeR::cpm	(e))
  rownames(cpmtable) = substr(rownames(cpmtable),1,15)
  return(cpmtable)
}

load_mutants = function(){
  file = mut_tpm_file
  #file = mut_fpkm_file
  raw = read.delim(file,row.names = 1)
  keep = apply(raw,1,function(x) sum(x) > 0.5)
  rpkmtable = raw[keep,]
  rownames(rpkmtable) = substr(rownames(rpkmtable),1,15)
  return(rpkmtable)  
}

load_mutcounts = function(){
  file = mut_count_file
  raw = read.delim(file,row.names = 1)
  keep = rowSums(edgeR::cpm(raw)>1) >= 2
  counttable = raw[keep, ]
  rownames(counttable) = substr(rownames(counttable),1,15)
  return(counttable)
}

load_mutcpm = function(){
  file = mut_count_file
  raw = read.delim(file,row.names = 1)
  keep = rowSums(edgeR::cpm(raw)>1) >= 2
  counttable = raw[keep, ]
  e = DGEList(counts=counttable)
  cpmtable = data.frame(cpm(e))
  rownames(cpmtable) = substr(rownames(cpmtable),1,15)
  return(cpmtable)
}

load_main = function(){
  file = mut_tpm_file
  #file = mut_fpkm_file
  raw = read.delim(file,row.names = 1)
  raw = select(raw,contains("es03"),contains("smad4"),contains("zic2"))
  keep = apply(raw,1,function(x) sum(x) > 0.5)
  rpkmtable = raw[keep,]
  rownames(rpkmtable) = substr(rownames(rpkmtable),1,15)
  return(rpkmtable)  
}

load_maincounts = function(){
  file = mut_count_file
  raw = read.delim(file,row.names = 1)
  raw = select(raw,contains("es03"),contains("smad4"),contains("zic2"))
  keep = rowSums(edgeR::cpm(raw)>1) >= 2
  counttable = raw[keep, ]
  rownames(counttable) = substr(rownames(counttable),1,15)
  return(counttable)
}

load_maincpm = function(){
  file = mut_count_file
  raw = read.delim(file,row.names = 1)
  raw = select(raw,contains("es03"),contains("smad4"),contains("zic2"))
  keep = rowSums(edgeR::cpm(raw)>1) >= 2
  counttable = raw[keep, ]
  e = DGEList(counts=counttable)
  cpmtable = data.frame(cpm(e))
  rownames(cpmtable) = substr(rownames(cpmtable),1,15)
  return(cpmtable)
}



load_smad4 = function(){
  file = mut_tpm_file
  #file = mut_fpkm_file
  raw = read.delim(file,row.names = 1)
  raw = select(raw,contains("es03"),contains("smad4"))
  keep = apply(raw,1,function(x) sum(x) > 0.5)
  rpkmtable = raw[keep,]
  rownames(rpkmtable) = substr(rownames(rpkmtable),1,15)
  return(rpkmtable)  
}

load_smad4counts = function(){
  file = mut_count_file
  raw = read.delim(file,row.names = 1)
  raw = select(raw,contains("es03"),contains("smad4"))
  keep = rowSums(edgeR::cpm(raw)>1) >= 2
  counttable = raw[keep, ]
  rownames(counttable) = substr(rownames(counttable),1,15)
  return(counttable)
}

load_smad4cpm = function(){
  file = mut_count_file
  raw = read.delim(file,row.names = 1)
  raw = select(raw,contains("es03"),contains("smad4"))
  keep = rowSums(edgeR::cpm(raw)>1) >= 2
  counttable = raw[keep, ]
  e = DGEList(counts=counttable)
  cpmtable = data.frame(edgeR::cpm(e))
  rownames(cpmtable) = substr(rownames(cpmtable),1,15)
  return(cpmtable)
}

load_zic2 = function(){
  file = mut_tpm_file
  #file = mut_fpkm_file
  raw = read.delim(file,row.names = 1)
  raw = select(raw,contains("es03"),contains("zic2"))
  keep = apply(raw,1,function(x) sum(x) > 0.5)
  rpkmtable = raw[keep,]
  rownames(rpkmtable) = substr(rownames(rpkmtable),1,15)
  return(rpkmtable)  
}

load_zic2counts = function(){
  file = mut_count_file
  raw = read.delim(file,row.names = 1)
  raw = select(raw,contains("es03"),contains("zic2"))
  keep = rowSums(edgeR::cpm(raw)>1) >= 2
  counttable = raw[keep, ]
  rownames(counttable) = substr(rownames(counttable),1,15)
  return(counttable)
}

load_zic2cpm = function(){
  file = mut_count_file
  raw = read.delim(file,row.names = 1)
  raw = select(raw,contains("es03"),contains("zic2"))
  keep = rowSums(edgeR::cpm(raw)>1) >= 2
  counttable = raw[keep, ]
  e = DGEList(counts=counttable)
  cpmtable = data.frame(edgeR::cpm(e))
  rownames(cpmtable) = substr(rownames(cpmtable),1,15)
  return(cpmtable)
}

load_day6 = function(){
  file = mut_tpm_file
  #file = mut_fpkm_file
  raw = read.delim(file,row.names = 1)
  raw = select(raw,contains(".d6."))
  keep = apply(raw,1,function(x) sum(x) > 0.5)
  rpkmtable = raw[keep,]
  rownames(rpkmtable) = substr(rownames(rpkmtable),1,15)
  return(rpkmtable)  
}

load_day6counts = function(){
  file = mut_count_file
  raw = read.delim(file,row.names = 1)
  raw = select(raw,contains(".d6."))
  keep = rowSums(edgeR::cpm(raw)>1) >= 2
  counttable = raw[keep, ]
  rownames(counttable) = substr(rownames(counttable),1,15)
  return(counttable)
}

load_day6cpm = function(){
  file = mut_count_file
  raw = read.delim(file,row.names = 1)
  raw = select(raw,contains(".d6."))
  keep = rowSums(edgeR::cpm(raw)>1) >= 2
  counttable = raw[keep, ]
  e = DGEList(counts=counttable)
  cpmtable = data.frame(edgeR::cpm(e))
  rownames(cpmtable) = substr(rownames(cpmtable),1,15)
  return(cpmtable)
}

load_msc_counts = function(){
  file = msc_count_file
  raw = read.delim(file,row.names = 1)
  keep = rowSums(edgeR::cpm(raw)>1) >= 2
  counttable = raw[keep, ]
  rownames(counttable) = substr(rownames(counttable),1,15)
  return(counttable)
}
#======================begin of data loading functions===================#
