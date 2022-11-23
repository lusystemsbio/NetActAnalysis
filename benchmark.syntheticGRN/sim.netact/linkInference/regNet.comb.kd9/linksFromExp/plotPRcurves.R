rm(list = ls())

inputdir <-  './data/'
regNetDir <- '../../regNet/data/' 

figdir <- "./figs/" 
dir.create(figdir) 

WIDTH <- 4
HEIGHT <- 8


# Precision recall using inferred links with NO DPI call
#------------------------------------------------------- 
ppvRcall <- read.table(file = paste(inputdir, 'precision-recall-noDPI.tsv', sep = ''), header = TRUE)

miTSH.list <- as.numeric(rownames(ppvRcall))

figname <- paste(figdir, 'performance-noDPI-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(3,1))
par(mar=c(2.5, 3.0, 1.5, 1.0))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space


# precision0-recall curve
plot(ppvRcall$rCall, ppvRcall$ppv, 
     lty=1, pch=19, 
     xlim = c(0, 0.3),
     ylim = c(0, 1), 
     xlab = 'recall', ylab = 'precision') 
lines(ppvRcall$rCall, ppvRcall$ppv, lty=2) 

performVal <- ppvRcall$rCall
plot(performVal,  
     lty=1, pch=19, 
     ylim = c(0, 1), xlab = 'mutual information threshold', ylab = 'recall', 
     xaxt = "n")  
axis(1,                         # Define x-axis manually
     at = 1:length(miTSH.list),
     labels = miTSH.list)
lines(performVal, lty=2) 

performVal <- ppvRcall$ppv
plot(performVal, #xlim = c(1, length(rCall.list)), 
     lty=1, pch=19, 
     ylim = c(0, 1), xlab = 'mutual information threshold', ylab = 'precision', 
     xaxt = "n")  
axis(1,                         # Define x-axis manually
     at = 1:length(miTSH.list),
     labels = miTSH.list)
lines(performVal, lty=2)  

dev.off() 




# Precision recall using inferred links with DPI call
#-----------------------------------------------------
ppvRcall <- read.table(file = paste(inputdir, 'precision-recall-DPI.tsv', sep = ''), header = TRUE)


XLIMIT <- c(0, max(miTSH.list))

figname <- paste(figdir, 'performance-DPI-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(3,1))
par(mar=c(2.5, 3.0, 1.5, 1.0))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space




# precision0-recall curve
plot(ppvRcall$rCall, ppvRcall$ppv, 
     lty=1, pch=19, 
     xlim = XLIMIT,
     ylim = c(0, 1), 
     xlab = 'recall', ylab = 'precision') 
lines(ppvRcall$rCall, ppvRcall$ppv, lty=2) 

performVal <- ppvRcall$rCall
plot(performVal,  
     lty=1, pch=19, 
     ylim = XLIMIT, 
     xlab = 'mutual information threshold', ylab = 'recall', 
     xaxt = "n")  
axis(1,                         # Define x-axis manually
     at = 1:length(miTSH.list),
     labels = miTSH.list)
lines(performVal, lty=2) 

performVal <- ppvRcall$ppv
plot(performVal, #xlim = c(1, length(rCall.list)), 
     lty=1, pch=19, 
     ylim = c(0, 1), xlab = 'mutual information threshold', ylab = 'precision', 
     xaxt = "n")  
axis(1,                         # Define x-axis manually
     at = 1:length(miTSH.list),
     labels = miTSH.list)
lines(performVal, lty=2)  

dev.off() 

