rm(list = ls()) 

#figdir <- "./figs/" 
figdir <- "./" 
#dir.create(figdir) 

ppvRcall.list <- list()

# NetAct - link inference is performed using NetAct activities  
# regulon db with the following perturbations: 0, 5, 10, 15 
ppvRcall.list[['ppvRcall.netact.regdb00']] <- read.table("../../../../sim.netact/linkInference/regNet.wt/linksFromNetAct/regdb00/data/prValues.tsv", sep = '\t')
ppvRcall.list[['ppvRcall.netact.regdb05']] <- read.table("../../../../sim.netact/linkInference/regNet.wt/linksFromNetAct/regdb05/data/prValues.tsv", sep = '\t')
ppvRcall.list[['ppvRcall.netact.regdb10']] <- read.table("../../../../sim.netact/linkInference/regNet.wt/linksFromNetAct/regdb10/data/prValues.tsv", sep = '\t')
ppvRcall.list[['ppvRcall.netact.regdb15']] <- read.table("../../../../sim.netact/linkInference/regNet.wt/linksFromNetAct/regdb15/data/prValues.tsv", sep = '\t')

# grnboost2.aucell - link inference is performed using AUCell activities
# regulon db with the following perturbations: 0, 5, 10, 15 
ppvRcall.list[['ppvRcall.grnboost2.aucell.regdb00']] <- read.table("../../../../sim.grnboost2/regNet.wt/linksFromAUcellAct/regdb00/data/precision-recall.tsv", sep = '\t')
ppvRcall.list[['ppvRcall.grnboost2.aucell.regdb05']] <- read.table("../../../../sim.grnboost2/regNet.wt/linksFromAUcellAct/regdb05/data/precision-recall.tsv", sep = '\t')
ppvRcall.list[['ppvRcall.grnboost2.aucell.regdb10']] <- read.table("../../../../sim.grnboost2/regNet.wt/linksFromAUcellAct/regdb10/data/precision-recall.tsv", sep = '\t')
ppvRcall.list[['ppvRcall.grnboost2.aucell.regdb15']] <- read.table("../../../../sim.grnboost2/regNet.wt/linksFromAUcellAct/regdb15/data/precision-recall.tsv", sep = '\t')


saveRDS(ppvRcall.list, file = paste0(figdir, 'ppvRcall.list.rds'))


length(ppvRcall.list)
method.names <- c('netact(regdb00)',  'netact(regdb05)', 'netact(regdb10)', 'netact(regdb15)', 
                  'grnboost2+aucell(regdb00)', 'grnboost2+aucell(regdb05)', 'grnboost2+aucell(regdb10)', 
                  'grnboost2+aucell(regdb15)'
                  
                  ) 

length(method.names)
names(method.names) <- names(ppvRcall.list)
names(ppvRcall.list)

# R colors: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf 
color.all <- c('black', 'red', 'blue', 'magenta', 'darkcyan', 'orange', 'brown', 'cyan','green')

length(color.all)
colorv <- color.all[1:length(ppvRcall.list)]
length(colorv)
names(colorv) <- names(ppvRcall.list)


WIDTH <- 6
HEIGHT <- 6
figname <- paste(figdir, 'prCurves-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(1,1)) 
par(mar=c(2.5, 3.0, 1.5, 1.0))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space

XLIMIT <- c(0, 1.0) #c(0, 0.10) #
YLIMIT <-  c(0, 1.0)

m.name <- names(ppvRcall.list)[1]
ppvRcall <- ppvRcall.list[[m.name]] 
colorCur <- colorv[m.name]
plot(ppvRcall[,'rCall'], ppvRcall[,'ppv'], 
     lty=1, pch=19,
     col=colorCur, 
     xlim = XLIMIT,
     ylim = YLIMIT,
     xlab = 'recall', ylab = 'precision', cex.lab=1.2)
lines(ppvRcall[,'rCall'], ppvRcall[,'ppv'], lty=2, col=colorCur)

# from second onward 
for(m.name in names(ppvRcall.list)[2:length(ppvRcall.list)]){
  print(m.name) 
  ppvRcall <-ppvRcall.list[[m.name]]
  colorCur <- colorv[m.name]
  points(ppvRcall[,'rCall'], ppvRcall[,'ppv'], 
         col=colorCur,
         lty=1, pch=19,
         xlim = c(0, 1),
         ylim = c(0, 1)) 
  lines(ppvRcall[,'rCall'], ppvRcall[,'ppv'], lty=2, col=colorCur)
  
  #break
}
#legend(x=0.45, y=0.45, legend = method.names, col = colorv, lwd = 2.5, bty = "n", cex=1.2)
dev.off() 

as.character(method.names) 
as.character(colorv)

figname <- paste(figdir, 'prCurve.legends-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(1,1)) 
par(mar=c(2.5, 3.0, 1.5, 1.0))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space

m.name <- names(ppvRcall.list)[1]
ppvRcall <- ppvRcall.list[[m.name]] 
colorCur <- colorv[m.name]
plot(ppvRcall[,'rCall'], ppvRcall[,'ppv'], 
     lty=1, pch=19, 
     col='white', 
     xlim = XLIMIT,
     ylim = YLIMIT,
     xlab = 'recall', ylab = 'precision')
legend(x=0.01, y=0.60, legend = method.names, col = colorv, lwd = 2.5, bty = "n")
dev.off() 
