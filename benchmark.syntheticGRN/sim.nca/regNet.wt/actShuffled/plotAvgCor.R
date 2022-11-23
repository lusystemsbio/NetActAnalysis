rm(list=ls())

data_dir <- figdir <- './cor/'

# Load correlations
#-------------------- 
cor.avg <- readRDS(file = paste(data_dir, './avg.cor.rds', sep = '')) 
cor.unperturbed <- read.csv('../actNoise.00/cor/avg.cor.nopert.csv', header = T)

# Plot correlations
#--------------------- 
#XLIMIT <- c(0.1, 0.45)
XLIMIT <- c(0.00, 1.00)
YLIMIT <- c(0, 300)

figname <- paste(figdir, 'hist.cor.avg.pdf', sep = '') 
pdf(file = figname, width=4, height=4, paper='special', onefile = TRUE)
par(mfrow=c(1,1)) 
par(mar=c(2.0, 2.0, 2.5, 0.5))  # bottom, left, top, right
hist(cor.avg,  #breaks = breaks, 
     #main = 'Avgerage correlations\ngene labels (racipe expr) shuffled', 
     main = '', 
     xlab = '', xlim = XLIMIT, ylim = YLIMIT)
abline(v=median(cor.avg), col='blue', lwd=2) 
abline(v=cor.unperturbed, col='red', lwd=2) 
dev.off()
