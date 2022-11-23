rm(list = ls()) 

#figdir <- "./figs/" 
figdir <- "./" 
#dir.create(figdir) 

ppvRcall.list <- list()

# NetAct - link inference is performed using NetAct activities  
# regulon db with the following perturbations: 5, 10, 15 i.e., 25%, 50%, 75% targets
ppvRcall.list[['ppvRcall.netact.regdb05']] <- read.table("../../../../sim.netact/linkInference/regNet.wt/linksFromNetAct/regdb05/data/stats4newLinks.tsv", sep = '\t')
ppvRcall.list[['ppvRcall.netact.regdb10']] <- read.table("../../../../sim.netact/linkInference/regNet.wt/linksFromNetAct/regdb10/data/stats4newLinks.tsv", sep = '\t')
ppvRcall.list[['ppvRcall.netact.regdb15']] <- read.table("../../../../sim.netact/linkInference/regNet.wt/linksFromNetAct/regdb15/data/stats4newLinks.tsv", sep = '\t')
saveRDS(ppvRcall.list, file = paste0(figdir, 'ppvRcall.list.newlinks.rds'))

length(ppvRcall.list)
method.names <- c('netact(regdb25)', 'netact(regdb50)', 'netact(regdb75)'#, 
                  #'genie3+aucell(regdb00)', 'genie3+aucell(regdb05)', 'genie3+aucell(regdb10)', 'genie3+aucell(regdb15)'
                  ) 

length(method.names)
names(method.names) <- names(ppvRcall.list)
names(ppvRcall.list)

# R colors: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf 
#color.all <- c('black', 'red', 'blue', 'magenta', 'darkcyan', 'orange', 'brown', 'cyan','green')

color.all <- c('orange', 'brown', 'cyan','green')

length(color.all)
colorv <- color.all[1:length(ppvRcall.list)]
length(colorv)
names(colorv) <- names(ppvRcall.list)

CEX.LAB <- 1.8
CEX.AXIS <- 1.2
#PCH.VAL <- 3
#PCH.VAL <- c(19, 10, 11, 12) 
PCH.VAL <- c(10, 11, 12) 
names(PCH.VAL) <- names(ppvRcall.list)
LTY.VAL <- 5
CEX.VAL <- 1.8

WIDTH <- 6
HEIGHT <- 6
figname <- paste(figdir, 'prCurves-newlinks-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(1,1)) 
par(mar=c(3.5, 3.5, 1.5, 1.0)+0.1)  # bottom, left, top, right
par(mgp=c(2.5, 0.8, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space

XLIMIT <- c(0, 1.0) #c(0, 0.10) #
YLIMIT <-  c(0, 1.0)

m.name <- names(ppvRcall.list)[1]
ppvRcall <- ppvRcall.list[[m.name]] 
colorCur <- colorv[m.name]
plot(ppvRcall[,'rCall'], ppvRcall[,'ppv'], 
     lty=1, pch=PCH.VAL[m.name],
     col=colorCur, 
     xlim = XLIMIT,
     ylim = YLIMIT,
     cex=CEX.VAL,
     xaxt='n',
     yaxt='n',
     xlab = 'recall', ylab = 'precision', cex.lab=CEX.LAB, cex.axis=CEX.AXIS)
lines(ppvRcall[,'rCall'], ppvRcall[,'ppv'], lty=LTY.VAL, col=colorCur)

# from second onward 
for(m.name in names(ppvRcall.list)[2:length(ppvRcall.list)]){
  print(m.name) 
  ppvRcall <-ppvRcall.list[[m.name]]
  colorCur <- colorv[m.name]
  points(ppvRcall[,'rCall'], ppvRcall[,'ppv'], 
         col=colorCur,
         lty=1, pch=PCH.VAL[m.name],
         cex=CEX.VAL,
         xlim = c(0, 1),
         ylim = c(0, 1), cex.lab=CEX.LAB, cex.axis=CEX.AXIS
         ) 
  lines(ppvRcall[,'rCall'], ppvRcall[,'ppv'], lty=LTY.VAL, col=colorCur)
  
  #break
}
axis(labels=NA,side=1,tck=-0.04,at=seq(0, 1, by=0.2))
axis(labels=NA,side=2,tck=-0.04,at=seq(0, 1, by=0.2)) 
#legend(x=0.45, y=0.45, legend = method.names, col = colorv, lwd = 2.5, bty = "n", cex=1.2)
dev.off() 

as.character(method.names) 
as.character(colorv)

figname <- paste(figdir, 'prCurve.legends-newlinks-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(1,1)) 
par(mar=c(2.5, 3.0, 1.5, 1.0))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space

m.name <- names(ppvRcall.list)[1]
ppvRcall <- ppvRcall.list[[m.name]] 
colorCur <- colorv[m.name]
plot(ppvRcall[,'rCall'], ppvRcall[,'ppv'], 
     lty=1, pch=PCH.VAL[m.name], 
     col='white', 
     xlim = XLIMIT,
     ylim = YLIMIT,
     xlab = 'recall', ylab = 'precision')
legend(x=0.01, y=0.60, legend = method.names, col = colorv, pch = PCH.VAL, bty = "n")
dev.off() 

