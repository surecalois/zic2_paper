source('/Users/jiexu/Desktop/test/my_functions.R')
source('/Users/jiexu/Desktop/test/zic2_functions.R')


counttable = load_zic2counts()
smad4_MDSedgeR(counttable)

raw = load_zic2counts()
fit.result = make.gmlFit(raw) 
desc = make.design(colnames(raw))
dz = c("(zic2c1.d6+zic2c2.d6)/2-es03.d6")
my.contrasts = makeContrasts(contrasts = dz, levels = desc)

dT.A = diff_table(fit.result,my.contrasts[,1])
A = diff_tops(dT.A,up = 100,down = 0)
A = A[order(A$FDR),]

lrt <- glmLRT(fit.result, contrast=my.contrasts[,1])
diff_mdplot(lrt)
  
#======================display a gene list=====================#
raw = load_zic2cpm()
#raw = load_zic2()
#genelist = load_list('CMD')
genelist = load_list()
expr.table = norm.by.max(raw[genelist,])
expr.table = prep_table(expr.table)
wrap_heatmap_list(expr.table,c('es03','zic2c1','zic2c2'))


#=====================kmeans============================#
raw = load_zic2cpm()
a = raw
a = a[filter.low(a,1) ,]
a = a[filter.similar(a, .9) ,]
a = a[filter.sd(a) ,]

b = kmeans_it(a,cluster_num = 12)
h = prep_table(b[,names(a)],check.complete = T)
groups = c('es03','zic2c1','zic2c2')
wrap_heatmap_long(h,groups = groups)
