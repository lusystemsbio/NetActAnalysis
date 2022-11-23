rm(list=ls())

NO.SIG.DIGITS <- 4

outdir <- './cor/'
dir.create(outdir)

# Load activities for current method
act.list <- readRDS(file = './act/act.all.rds') 
act.curMethod <- act.list[[1]]

# Load racipe activities 
act.racipe.df <- read.table(file = '../../../sim.racipe/sim.racipe.wt/net30tf.states.txt', header = TRUE) 
rownames(act.racipe.df) <- paste('M', act.racipe.df$MODEL_NO, sep = '')
act.racipe <- t(act.racipe.df[, rownames(act.curMethod)]) # retain racipe act for TFs only


# Calculate correlations
cor.avg.v <- sapply(1:length(act.list), function(item_no)
  mean(abs(diag(cor(t(act.racipe), t(act.list[[item_no]]), method = 'spearman'))))
)
 
saveRDS(cor.avg.v, paste(outdir, 'cor.avg.rds', sep = ''))
