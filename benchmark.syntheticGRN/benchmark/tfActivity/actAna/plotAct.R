rm(list = ls())

figdir <- './figs/'
dir.create(figdir)

regdb <- readRDS(file = '/Users/a.katebi/research/netact/benchmark/regDBs/regDB.nopert.rds')
targets <- unique(unlist(regdb))
tfs <- names(regdb)
length(targets)

act.racipe.df <- read.table(file = '/Users/a.katebi/research/netact/benchmark/sim.racipe/sim.racipe.wt/net30tf.states.txt', header = TRUE)
# act.racipe.df <- act.racipe.df[, tfs]
rownames(act.racipe.df) <- paste('M', act.racipe.df$MODEL_NO, sep = '')
act.racipe.df$MODEL_NO <- NULL; act.racipe.df$NO_STATES <- NULL; act.racipe.df$STATE_NO <- NULL
act.racipe.df <- t(act.racipe.df)

act.netact.df <- read.csv(file = '/Users/a.katebi/research/netact/benchmark/sim.netact/regNet.wt/actNoise.00/act/act.regDB.0.csv', header = TRUE, row.names = 1)
act.netact.df <- as.matrix(act.netact.df)

act.nca.df <- read.csv(file = '/Users/a.katebi/research/netact/benchmark/sim.nca/regNet.wt/actNoise.00/act/act.regDB.0.csv', header = TRUE, row.names = 1)
act.nca.df <- as.matrix(act.nca.df)

act.viper.df <- read.csv(file = '/Users/a.katebi/research/netact/benchmark/sim.viper/regNet.wt/actNoise.00/act/act.regDB.0.csv', header = TRUE, row.names = 1)
act.viper.df <- as.matrix(act.viper.df)

act.aucell.1.df <- read.csv(file = '/Users/a.katebi/research/netact/benchmark/sim.aucell.1/regNet.wt/actNoise.00/act/act.regDB.0.csv', header = TRUE, row.names = 1)
act.aucell.2.df <- read.csv(file = '/Users/a.katebi/research/netact/benchmark/sim.aucell.2/regNet.wt/actNoise.00/act/act.regDB.0.csv', header = TRUE, row.names = 1)
act.aucell.1.df <- as.matrix(act.aucell.1.df)
act.aucell.2.df <- as.matrix(act.aucell.2.df)


scale_data <- function(x) y <- (x-min(x))/(max(x)-min(x)) # bring data between 0 and 1

act.racipe.df <- scale_data(act.racipe.df)
act.nca.df <- scale_data(act.nca.df)
act.netact.df <- scale_data(act.netact.df)
act.viper.df <- scale_data(act.viper.df)
act.aucell.1.df <- scale_data(act.aucell.1.df)
act.aucell.2.df <- scale_data(act.aucell.2.df)

YLIMIT <- c(0.0, 1.1) 
XLIMIT <- c(1, 83) 
PCH_VAL <- 20
colorv <- c('red', 'blue', 'cyan', 'orange', 'forestgreen')
legend_names <- c('RACIPE', 'NetAct', 'VIPER', 'AUCell 1', 'AUCell 2')


WIDTH <- 8
HEIGHT <- 4
figname <- paste(figdir, 'ActByMethod-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)

par(mfrow=c(1,1))
par(mar=c(3.5, 3.5, 1.5, 1.0))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(1,1,1,1)) # all sides have 3 lines of space

for(tf in tfs[1:9]){
  print(tf)  
  plot(1:length(XLIMIT), xlim=XLIMIT, 
       ylim=YLIMIT, main = tf,
       pch='', #xaxt='n', 
       #lwd=0.01,
       cex=0.01, xlab = 'Model No', ylab = 'Activity')
  points(act.racipe.df[tf,], col=colorv[1], pch = PCH_VAL)
  lines(act.racipe.df[tf,], col=colorv[1], lty=3)
  
  points(act.netact.df[tf,], col=colorv[2], pch = PCH_VAL)
  lines(act.netact.df[tf,], col=colorv[2], lty=3)
  
  points(act.viper.df[tf,], col=colorv[3], pch = PCH_VAL)
  lines(act.viper.df[tf,], col=colorv[3], lty=3)
  
  points(act.aucell.1.df[tf,], col=colorv[4], pch = PCH_VAL)
  lines(act.aucell.1.df[tf,], col=colorv[4], lty=3)
  
  points(act.aucell.2.df[tf,], col=colorv[5], pch = PCH_VAL)
  lines(act.aucell.2.df[tf,], col=colorv[5], lty=3)
  
  abline(v=1:ncol(act.racipe.df), lwd=0.15)
  legend(1.2, 1.15, legend = legend_names, col=colorv, lty=3, pch = PCH_VAL, horiz = TRUE, bty = 'n')
  #break()
} 
dev.off()

# scatter plot for TF activities between methods
WIDTH <- 4
HEIGHT <- 6
figname <- paste(figdir, 'ActScatter-RACIPE-vs-otherMethods-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(4,1))
par(mar=c(0.1, 4.0, 1.5, 1.5))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(4,1,1,1)) # all sides have 3 lines of space
for(tf in tfs[1:9]){ 
  print(tf)
  plot(act.racipe.df[tf,], act.netact.df[tf,], xlab = '', ylab = 'Netact', xaxt='n')
  title(main = tf)
  plot(act.racipe.df[tf,], act.viper.df[tf,], xlab = '', ylab = 'VIPER', xaxt='n')
  plot(act.racipe.df[tf,], act.aucell.1.df[tf,], xlab = '', ylab = 'AUCell.1', xaxt='n')
  plot(act.racipe.df[tf,], act.aucell.2.df[tf,], xlab = 'RACIPE', ylab = 'AUCell.2')
  
  #title(xlab = 'RACIPE',main = 'RACIPE vs other methods')
}
dev.off()


cordf <- as.data.frame(matrix(nrow = length(tfs), ncol = 4))
rownames(cordf) <- tfs
colnames(cordf) <- c('netact' , 'viper', 'aucell.1', 'aucell.2')
for(tf in tfs){ 
  cordf[tf, 'netact'] <- cor(act.racipe.df[tf,], act.netact.df[tf,], method = 'spearman')
  cordf[tf, 'viper'] <- cor(act.racipe.df[tf,], act.viper.df[tf,], method = 'spearman')
  cordf[tf, 'aucell.1'] <- cor(act.racipe.df[tf,], act.aucell.1.df[tf,], method = 'spearman')
  cordf[tf, 'aucell.2'] <- cor(act.racipe.df[tf,], act.aucell.2.df[tf,], method = 'spearman')
}

tfOrder <- paste('tf', 1:30, sep = '')
tfOrder
cordf <- cordf[tfOrder,]
tfOrderLabel <- paste('TF', 1:30, sep = ' ')
methodNames <- c('Netact', 'VIPER', 'AUCell 1', 'AUCell 2')

WIDTH <- 8
HEIGHT <- 4
figname <- paste(figdir, 'Cor-RACIPE-vs-otherMethods-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(1,1))
par(mar=c(2.5, 4.0, 1.5, 1.5))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space
plot(1:length(tfs), xlim=c(1, 30), ylim=c(-1, 1.2), main = '', pch='', 
     xaxt='n', #lwd=0.01,
     cex=0.01, xlab = 'TF No', ylab = 'Cor')
points(cordf$netact, col=colorv[1], pch = PCH_VAL)
lines(cordf$netact, col=colorv[1], lty=3)

points(cordf$viper, col=colorv[2], pch = PCH_VAL)
lines(cordf$viper, col=colorv[2], lty=3) 

points(cordf$aucell.1, col=colorv[3], pch = PCH_VAL)
lines(cordf$aucell.1, col=colorv[3], lty=3)

points(cordf$aucell.2, col=colorv[4], pch = PCH_VAL)
lines(cordf$aucell.2, col=colorv[4], lty=3)

abline(v=1:length(tfs), lwd=0.1, lty=3)
axis(1, at=1:length(tfOrderLabel), las=1, cex.axis = 1.0, tck = -0.02) 
legend(1.2, 1.3, legend = methodNames, col=colorv, lty=3, pch = PCH_VAL, horiz = TRUE, bty = 'n')
dev.off()
