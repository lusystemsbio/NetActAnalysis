rm(list=ls()) 

data_dir <- figdir <- './cor/'

cor.avg <- readRDS(paste(data_dir, './cor.avg.rds', sep=''))
cor.unperturbed <- unlist(read.csv('../actNoise.00/cor/avg.cor.nopert.csv', header = T))

# Plot correlations
#--------------------- 
XLIMIT <- c(0.0, 0.5) 
YLIMIT <- c(0, 300)
#breaks = seq(0.5, 1, by=0.1)

figname <- paste(figdir, 'hist.avg.cor.pdf', sep = '') 
pdf(file = figname, width=4, height=4, paper='special', onefile = TRUE)
par(mfrow=c(1,1))
par(mar=c(2.0, 2.0, 1.5, 0.5))  # bottom, left, top, right

hist(cor.avg,  #breaks = breaks, 
     #main = 'Avg correlations: gene labels (racipe expressions) shuffled', 
     main = '', 
     xlab = '', xlim = XLIMIT, ylim = YLIMIT)  
abline(v=cor.unperturbed, col='red', lwd=2)
abline(v=mean(cor.avg), col='blue', lwd=2)
dev.off()
