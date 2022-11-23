rm(list=ls())

# Load correlations
#-------------------
data_dir <- './cor/'

fname.cor <- paste(data_dir, './avg.cor.pert.rds', sep='')
cor.avg.list <- readRDS(fname.cor)

cor.unperturbed <- read.csv(paste(data_dir, "avg.cor.nopert.csv", sep = ''), header = T ) 
cor.unperturbed <- unlist(cor.unperturbed)

cor.avg.regDB.5 <- mean(unlist(cor.avg.list$regdb.5))
print(cor.avg.regDB.5)

cor.avg.regDB.10 <- mean(unlist(cor.avg.list$regdb.10))
print(cor.avg.regDB.10)

cor.avg.regDB.15 <- mean(unlist(cor.avg.list$regdb.15))
print(cor.avg.regDB.15)

# Plot correlations
#--------------------- 
XLIMIT <- c(0, 1.0)
#XLIMIT <- c(0.5, 0.90)
#XLIMIT <- c(0.2, 0.50) 

#XLIMIT <- c(0.2, 0.55)
#YLIMIT <- c(0, 25) 
YLIMIT <- c(0, 30) 
#breaks = seq(0.5, 1, by=0.1) 

#figdir <- paste(data_dir, "figs/", sep='') 
#dir.create(figdir) 
figdir <- data_dir

figname <- paste(figdir, 'dist.cor.avg.pdf', sep = '') 
pdf(file = figname, width=6, height=6, 
    paper='special', onefile = TRUE)
par(mfrow=c(3,1))
par(mar=c(2.0, 2.0, 1.5, 0.5))  # bottom, left, top, right

hist(unlist(cor.avg.list$regdb.5),  #breaks = breaks, 
     main = 'Avgerage correlations: regDB 5',  
     xlim = XLIMIT, ylim = YLIMIT, 
     xlab = '')  
abline(v=cor.unperturbed, col='red', lwd=2)
abline(v=cor.avg.regDB.5, col='blue', lwd=2)


hist(unlist(cor.avg.list$regdb.10), #breaks = breaks, 
     main = 'regDB 10', 
     xlim = XLIMIT, ylim = YLIMIT, 
     xlab = '')  
abline(v=cor.unperturbed, col='red', lwd=2)
abline(v=cor.avg.regDB.10, col='blue', lwd=2)

hist(unlist(cor.avg.list$regdb.15), #breaks = breaks, 
     main = 'regDB 15',  
     xlim = XLIMIT, ylim = YLIMIT, 
     xlab = '')  
abline(v=cor.unperturbed, col='red', lwd=2)
abline(v=cor.avg.regDB.15, col='blue', lwd=2)
dev.off() 
