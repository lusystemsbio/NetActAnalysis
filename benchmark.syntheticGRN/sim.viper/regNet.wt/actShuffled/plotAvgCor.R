rm(list=ls())

data_dir <- './cor/'
figdir <- data_dir


# Load correlations
#--------------------
cor.avg <- readRDS(paste(data_dir, './cor.avg.rds', sep='')) # avg cor from randomized racipe expressions

cor.unperturbed <- read.csv('../actNoise.00/cor/avg.cor.nopert.csv', header = T) # avg cor from WT racipe expressions
cor.unperturbed <- unlist(cor.unperturbed)

# Plot correlations
#------------------- 
#XLIMIT <- c(0, 0.40)
XLIMIT <- c(0.00, 1.00)
YLIMIT <- c(0, 1700)
#breaks = seq(0.5, 1, by=0.1)

pdf(file = paste(figdir, 'dist.cor.avg.pdf', sep = '') , width=4, height=4, 
    paper='special', onefile = TRUE)
par(mfrow=c(1,1))
par(mar=c(2.0, 2.0, 1.5, 0.5))  # bottom, left, top, right

hist(cor.avg,  #breaks = breaks, 
     #main = 'Avg correlations: gene labels shuffled', xlab = '', 
     main = '', xlab = '', 
     xlim = XLIMIT, ylim = YLIMIT)  
abline(v=cor.unperturbed, col='red', lwd=2)
abline(v=mean(cor.avg), col='blue', lwd=2)
dev.off() 
