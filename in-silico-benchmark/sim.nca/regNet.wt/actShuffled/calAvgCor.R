rm(list=ls())

NO.SIG.DIGITS <- 4

outdir <- './cor/'
dir.create(outdir)

act.all <- readRDS(file = './act/act.all.rds') 

# Load racipe activities 
act.racipe.df <- read.table(file = '../../../sim.racipe/sim.racipe.wt/net30tf.states.txt', header = TRUE) 
rownames(act.racipe.df) <- paste('M', act.racipe.df$MODEL_NO, sep = '') 
act.racipe <- t(act.racipe.df[, rownames( act.all[[1]])]) # retain only TFs
sum(rownames(act.racipe) == rownames(act.all[[1]])) # check TF order

# First, calculate cor between racipe and nca activities for each TF and then calculate their mean
cor.avg <- sapply(1:length(act.all), function(k) mean(abs(diag(cor(t(act.racipe), t(act.all[[k]]))))))  
length(cor.avg)
saveRDS(cor.avg, file = paste(outdir, 'avg.cor.rds', sep = ''))
