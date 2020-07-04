#library(DESeq2)
library(edgeR)
library(ggplot2)

datafile = "/Users/jiexu/kclab/crispr_screen/screen/V3.5/res_v3.5/res1000-zlib2all.txt"
fout = F
#DESeq_outfile = "/Users/jiexu/kclab/crispr_screen/screen/V3.5/res_v3.5/mesp1hESC_DESeq_output.txt"
edgeR_outfile = "/Users/jiexu/kclab/crispr_screen/screen/V3.5/res_v3.5/mesp1hESC_edgeR_output.txt"
#merge_outfile = "/Users/jiexu/kclab/crispr_screen/screen/V3.5/res_v3.5/mesp1hESC_merge_output.txt"

raw <- read.table(datafile, header=F);
gp_rp=c("neg.A","pos.A","neg.B","pos.B","neg.C","pos.C")
gps=c("neg","pos")
rpn=c("A","B","C")

isl1.lib=paste("isl1",gp_rp,sep=".")
mesp1.lib=paste("mesp1",gp_rp,sep=".")
hesc.lib=paste("zlib2",rpn,sep=".")

names(raw)=c("chr","st","ed","sgRNA","gene","sig",isl1.lib,mesp1.lib,hesc.lib,"total","pam","seq")
rownames(raw)=raw[,"sgRNA"]

counttable = raw[,isl1.lib]
condition = c("neg","pos","neg","pos","neg","pos")

isl1.sum = data.frame(isl1.A = raw$isl1.neg.A+raw$isl1.pos.A,
                      isl1.B = raw$isl1.neg.B+raw$isl1.pos.B,
                      isl1.C = raw$isl1.neg.C+raw$isl1.pos.C)

mesp1.sum = data.frame(mesp1.A = raw$mesp1.neg.A+raw$mesp1.pos.A,
                       mesp1.B = raw$mesp1.neg.B+raw$mesp1.pos.B,
                       mesp1.C = raw$mesp1.neg.C+raw$mesp1.pos.C)

hESC.sum = data.frame(zlib2.A = 0.6*raw$zlib2.A+0.2*raw$zlib2.B+0.2*raw$zlib2.C,
                      zlib2.B = 0.2*raw$zlib2.A+0.6*raw$zlib2.B+0.2*raw$zlib2.C,
                      zlib2.C = 0.2*raw$zlib2.A+0.2*raw$zlib2.B+0.6*raw$zlib2.C)

counttable = cbind(mesp1.sum,isl1.sum)#round(hESC.sum))
rownames(counttable)=raw[,"sgRNA"]
condition = c("neg","neg","neg","pos","pos","pos")
#condition = c("mesp1","mesp1","mesp1","isl1","isl1","isl1")


keep = rowSums(counttable) >= 100
counttable = counttable[keep,]

meta = data.frame(
  row.names = colnames(counttable), 
  condition = condition,
  libType=c("single", "single","single", "single","single", "single")
)
meta$condition = relevel(meta$condition,ref = "pos");

###DESeq_
# dds = DESeqDataSetFromMatrix(countData = counttable,colData = meta,design = ~condition)
# dds = DESeq(dds)
# res = results(dds)
# DESeq_res=res
# res = cbind(counttable,res)
# resOrdered = res[order(res$padj),]
# DESeq_out = cbind("gene"=raw[rownames(resOrdered),"gene"],resOrdered)
# if(fout) write.table(as.data.frame(DESeq_out),sep="\t",quote=F,file=DESeq_outfile,row.names = F)

###edgeR_
condition <- relevel(factor(meta$condition), ref="pos")
libType <- factor(meta$libType)
e = DGEList(counts=counttable,group = condition)
e = estimateCommonDisp(e)
e = estimateTagwiseDisp(e)
et = exactTest(e)
res <- topTags(et, n=nrow(e),sort.by="none")$table
edgeR_res=res
res = cbind(counttable,res)
resOrdered = res[order(res$PValue),] #res[order(res$FDR),]
edgeR_out = cbind("gene"=raw[rownames(resOrdered),"gene"],resOrdered)

if(fout) write.table(as.data.frame(edgeR_out),sep="\t",quote=F,file=edgeR_outfile,row.names = F)

dt = decideTestsDGE(et,lfc=3,p.value = 1e-3)
detags = rownames(et$table)[as.logical(dt)]
des = et$table[detags,]
pbasic = ggplot()+ggtitle("mesp1 vs. isl1")#+ylim(c(-12,12))
p0 = geom_point(data = et$table,aes(y = logFC,x = logCPM),alpha = I(3/7),size = 2)
pd = geom_point(data = des,aes(y = logFC,x = logCPM),color = "red" ,size = 2)
print(pbasic+pd+p0)




###all
# all=cbind(counttable,edgeR_res,DESeq_res)
# all$sort = sqrt(all$padj*all$FDR)
# resOrdered = all[order(all$sort),]
# merge_out = cbind("gene"=raw[rownames(resOrdered),"gene"],resOrdered)
# if(fout) write.table(as.data.frame(merge_out),sep="\t",quote=F,file=merge_outfile,row.names = F)


#=======some MDS plot==============#
# test = raw[,c(isl1.lib,mesp1.lib,hesc.lib)]
# test_c = c("isl1.neg","isl1.pos","isl1.neg","isl1.pos","isl1.neg","isl1.pos",
#            "mesp1.neg","mesp1.pos","mesp1.neg","mesp1.pos","mesp1.neg","mesp1.pos",
#            "hESC","hESC","hESC")
# rep_c = c("A","A","B","B","C","C","A","A","B","B","C","C","A","B","C")
# keep = apply(test[,hesc.lib],1,function(x) min(x) > 100)
# #keep = apply(test,1,function(x) max(x) > 500)
# test = test[keep,]
# te = DGEList(counts=test,group = test_c)
# 
# z = plotMDS(te,method = "logFC",top = 5000)#"bcv""logFC"
# yp = data.frame(x = z$x,y = -z$y,type = test_c,rep = rep_c)
# yp$label = paste(yp$type,yp$rep,sep = ".")
# gg_base = ggplot(data = yp) + xlim(min(z$x)-1,max(z$x)+1)
# points = geom_point(aes(x = x, y = y,shape = rep,color = type),size = 5)
# points2 = geom_point(aes(x = x, y = y,shape = rep),color = "black",size = 6)
# label = geom_text(aes(x = x+.2, y = y,label = label),hjust = 0, nudge_x = 5,size = 3)
# print(gg_base+points2+points+label)

#keep = apply(test,1,function(x) min(x) > 0)
# sel = 10:500
# cmd.data = cpm(test[sel,],log = T)
# zzz = cmdscale(dist(cmd.data))
# zd = data.frame(x = zzz[,1],y = zzz[,2],gene = raw[sel,]$gene)
# zp = ggplot(data = zd,aes(x = x, y= y))+geom_point(alpha = I(1/5))+geom_text(aes(x = x,y = y, label = gene),size = 2)
#print(zp)

# file = "/Users/jiejiaxu/genome/annotation/list/focus_list.txt"
# listfile <- read.table(file, header=F);
# inlist = listfile$V1
# zdk = zd[zd$gene %in% inlist,]
# zpk = geom_text(data = zdk,aes(x = x,y = y, label = gene),color = "red",size = 2)
#print(zp+zpk)



