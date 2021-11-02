files = list.files('/project2/mstephens/dongyue/gtex/Yang/CoverageCounts')
f = files[3]
fs = strsplit(f,split='_')
gene_name = fs[[1]][1]
fs = strsplit(fs[[1]][2],split = ':')
chr = fs[[1]][1]

datax = read.table(paste('/project2/mstephens/dongyue/gtex/Yang/CoverageCounts/',f,sep=''),
                   header = TRUE)

idx = 4
y = (datax[,idx])
#y = rowSums(datax[,-c(1,2,3)])
indi_name = colnames(datax)[idx]

plot(y,col='grey50',type='p',pch=20,cex=0.3,ylab='counts',xlab='Chr11',xaxt = "n",main=paste(gene_name,indi_name,sep=', '))
xlab.idx = round(seq(1,nrow(datax),length.out = 5))
axis(1, at = xlab.idx,labels = datax[xlab.idx,2])

# fit smash poiss
fit.smash = smash.poiss(y)
plot(y,col='grey40',type='p',pch=20,cex=0.2,ylab='counts',xlab='Chr11',xaxt = "n",main=paste(gene_name,indi_name,sep=', '))
xlab.idx = round(seq(1,nrow(datax),length.out = 5))
axis(1, at = xlab.idx,labels = datax[xlab.idx,2])
lines(fit.smash,col='grey10')

fit.smash.gen = smash.gen.poiss(y,transformation = 'lik_expansion',
                                smooth_method = 'smash',smash.pm = F)
plot(y,col='grey50',type='p',pch=20,ylab='counts',
     xlab='Chr11',xaxt = "n",main=paste(gene_name,indi_name,sep=', '),
     cex=0.2)
xlab.idx = round(seq(1,nrow(datax),length.out = 5))
axis(1, at = xlab.idx,labels = datax[xlab.idx,2])
lines(fit.smash.gen$lambda.est,lwd=2,col=4)
