rm(list=ls())

outdir <- './cor/'
dir.create(outdir)

# Load calculated activities 
act.all <- readRDS(file = './act/act.all.rds') 

# Load racipe activities 
act.racipe.df <- read.table(file = '../../../sim.racipe/sim.racipe.wt/net30tf.states.txt', header = TRUE) 
act.racipe <- t(act.racipe.df[, rownames(act.all[[1]])]) # retain only TF activities

dim(act.racipe)
sum(rownames(act.racipe) == rownames(act.all[[1]])) # Check TF order

# Calculate correlations
cor.avg <- sapply(1:length(act.all), function(k) mean(abs(diag(cor(t(act.racipe), t(act.all[[k]]))))) )
length(cor.avg)
saveRDS(cor.avg, file = paste(outdir, 'cor.avg.rds', sep = ''))
