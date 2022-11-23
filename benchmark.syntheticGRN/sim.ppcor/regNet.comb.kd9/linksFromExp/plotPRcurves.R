rm(list = ls())

pcorTSH.list <- c(0.10, 0.15, 0.20, 0.25, 0.30)
miTSH.list <- pcorTSH.list

inputdir <-  './data/'
regNetDir <- '../../regNet/data/' 

figdir <- "./figs/" 
dir.create(figdir) 

WIDTH <- 4
HEIGHT <- 8

# Precision recall using inferred links with NO DPI call
#------------------------------------------------------- 
ppvRcall <- read.table(file = paste(inputdir, 'precision-recall.tsv', sep = ''), header = TRUE)


figname <- paste(figdir, 'performance-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(3,1))
par(mar=c(2.5, 3.0, 1.5, 1.0))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space

YLIMIT <- c(0, 1.0)

# precision0-recall curve
plot(ppvRcall$rCall, ppvRcall$ppv, 
     lty=1, pch=19, 
     xlim = c(0, 1.0),
     ylim = YLIMIT , 
     xlab = 'recall', ylab = 'precision') 
lines(ppvRcall$rCall, ppvRcall$ppv, lty=2) 

performVal <- ppvRcall$rCall
plot(performVal,  
     lty=1, pch=19, 
     ylim = YLIMIT, xlab = 'Partial corr threshold', ylab = 'recall', 
     xaxt = "n")  
axis(1,                         # Define x-axis manually
     at = 1:length(miTSH.list),
     labels = miTSH.list)
lines(performVal, lty=2) 

performVal <- ppvRcall$ppv
plot(performVal, #xlim = c(1, length(rCall.list)), 
     lty=1, pch=19, 
     ylim = YLIMIT , xlab = 'Partial corr threshold', ylab = 'precision', 
     xaxt = "n")  
axis(1,                         # Define x-axis manually
     at = 1:length(miTSH.list),
     labels = miTSH.list)
lines(performVal, lty=2)  

dev.off() 

