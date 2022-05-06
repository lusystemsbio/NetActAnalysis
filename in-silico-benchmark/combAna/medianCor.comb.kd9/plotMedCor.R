rm(list=ls())

figdir <- './figs/'
dir.create(figdir)

# Noise Level 00: load avg cor (no perturbation) and median cor (for three perturbed regulon)
avg_and_medianCor.list <- list() 

# Netact
methodDir <- '../../sim.Netact/regNet.comb.kd9/'
method_name <- strsplit(methodDir, split = 'sim.')
method_name <- method_name[[1]][2]
method_name <- strsplit(method_name, split = '/regNet.comb.kd9/')
method_name <- method_name[[1]]
method_name 

avgCor_noPert <- unlist(read.csv(file = paste0(methodDir, 'actNoise.00/cor/avg.cor.nopert.csv')))
meadianCor <- read.csv(file =  paste0(methodDir, 'actNoise.00/cor/median.cor.pert.csv')) 
avg_and_medianCor00 <- c(avgCor_noPert, meadianCor$median_cor)
names(avg_and_medianCor00) <- c('regdb.0', meadianCor$pert_level) 
avg_and_medianCor00
avg_and_medianCor.list[[method_name]] <- avg_and_medianCor00


# NCA
methodDir <- '../../sim.nca/regNet.comb.kd9/'
method_name <- strsplit(methodDir, split = 'sim.')
method_name <- method_name[[1]][2]
method_name <- strsplit(method_name, split = '/regNet.comb.kd9/')
method_name <- method_name[[1]]

avgCor_noPert <- unlist(read.csv(file = paste0(methodDir, 'actNoise.00/cor/avg.cor.nopert.csv')))
meadianCor <- read.csv(file =  paste0(methodDir, 'actNoise.00/cor/median.cor.pert.csv')) 
avg_and_medianCor00 <- c(avgCor_noPert, meadianCor$median_cor)
names(avg_and_medianCor00) <- c('regdb.0', meadianCor$pert_level) 
avg_and_medianCor00
avg_and_medianCor.list[[method_name]] <- avg_and_medianCor00


# VIPER
methodDir <- '../../sim.viper/regNet.comb.kd9/'
method_name <- strsplit(methodDir, split = 'sim.')
method_name <- method_name[[1]][2]
method_name <- strsplit(method_name, split = '/regNet.comb.kd9/')
method_name <- method_name[[1]]

avgCor_noPert <- unlist(read.csv(file = paste0(methodDir, 'actNoise.00/cor/avg.cor.nopert.csv')))
meadianCor <- read.csv(file =  paste0(methodDir, 'actNoise.00/cor/median.cor.pert.csv')) 
avg_and_medianCor00 <- c(avgCor_noPert, meadianCor$median_cor)
names(avg_and_medianCor00) <- c('regdb.0', meadianCor$pert_level) 
avg_and_medianCor00
avg_and_medianCor.list[[method_name]] <- avg_and_medianCor00


# AUCell 1: regulon DB modified for AUCell 
methodDir <- '../../sim.aucell.1/regNet.comb.kd9/'
method_name <- strsplit(methodDir, split = 'sim.')
method_name <- method_name[[1]][2]
method_name <- strsplit(method_name, split = '/regNet.comb.kd9/')
method_name <- method_name[[1]]

avgCor_noPert <- unlist(read.csv(file = paste0(methodDir, 'actNoise.00/cor/avg.cor.nopert.csv')))
meadianCor <- read.csv(file =  paste0(methodDir, 'actNoise.00/cor/median.cor.pert.csv')) 
avg_and_medianCor00 <- c(avgCor_noPert, meadianCor$median_cor)
names(avg_and_medianCor00) <- c('regdb.0', meadianCor$pert_level) 
avg_and_medianCor00
avg_and_medianCor.list[[method_name]] <- avg_and_medianCor00


# AUCell 2: base regulon DB   
methodDir <- '../../sim.aucell.2/regNet.comb.kd9/'
method_name <- strsplit(methodDir, split = 'sim.')
method_name <- method_name[[1]][2]
method_name <- strsplit(method_name, split = '/regNet.comb.kd9/')
method_name <- method_name[[1]]

avgCor_noPert <- unlist(read.csv(file = paste0(methodDir, 'actNoise.00/cor/avg.cor.nopert.csv')))
meadianCor <- read.csv(file =  paste0(methodDir, 'actNoise.00/cor/median.cor.pert.csv')) 
avg_and_medianCor00 <- c(avgCor_noPert, meadianCor$median_cor)
names(avg_and_medianCor00) <- c('regdb.0', meadianCor$pert_level) 
avg_and_medianCor00
avg_and_medianCor.list[[method_name]] <- avg_and_medianCor00


XLIMIT <- c(1, length(avg_and_medianCor.list$Netact))
YLIMIT <- c(0.0, 1.0) 
legend_names <- names(avg_and_medianCor.list)
pert_levels <- c('PERT 0', 'PERT 5', 'PERT 10', 'PERT 15')
#colorv <- c('red', 'blue', 'green', 'magenta')
#colorv <- c('red', 'blue', 'green', 'cyan') 
#colorv <- c('red', 'blue', 'cyan', 'maroon') 
colorv <- c('red', 'blue', 'cyan', 'orange', 'forestgreen') 

HEIGHT <- 4
WIDTH <- 4
figname <- paste(figdir, 'medianCor_allMethods','-', WIDTH, 'x', HEIGHT, '.pdf', sep = '') 
pdf(file = figname, width=WIDTH, height=HEIGHT, paper='special', onefile = TRUE)

par(mfrow=c(1,1)) 
par(mar=c(5.0, 4.5, 2.5, 1.0))  # bottom, left, top, right
plot(1:length(XLIMIT),xlim=XLIMIT, ylim=YLIMIT, pch='', xaxt='n', xlab = '', ylab = '')
axis(1, at=1:length(pert_levels), labels=pert_levels, las=2)
points(avg_and_medianCor.list$Netact, col=colorv[1], pch = 19)
lines(avg_and_medianCor.list$Netact, col=colorv[1], lty=3)

points(avg_and_medianCor.list$nca, col=colorv[2], pch = 19)
lines(avg_and_medianCor.list$nca, col=colorv[2], lty=3)

points(avg_and_medianCor.list$viper, col=colorv[3], pch = 19) 
lines(avg_and_medianCor.list$viper, col=colorv[3], lty=3)

points(avg_and_medianCor.list$aucell.1, col=colorv[4], pch = 19)
lines(avg_and_medianCor.list$aucell.1, col=colorv[4], lty=3)

points(avg_and_medianCor.list$aucell.2, col=colorv[5], pch = 19)
lines(avg_and_medianCor.list$aucell.2, col=colorv[5], lty=3)

#legend(1.2, 0.80, legend = legend_names, col=colorv, lty=3, pch = 19)
#legend(1.0, 0.78, legend = legend_names, col=colorv, lty=3, pch = 19, bty = 'n', ncol = 2)

#title(xlab = 'perturbation level', ylab = 'median correlation')
#title(ylab = 'mean/median cor')
#title(ylab = 'median accuracy')
title(ylab = 'median correlation')
dev.off()  
