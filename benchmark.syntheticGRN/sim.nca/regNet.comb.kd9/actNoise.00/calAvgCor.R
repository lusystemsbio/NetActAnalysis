rm(list=ls())

NO.SIG.DIGITS <- 4

outdir <- './cor/'
dir.create(outdir)

# Calculate cor for unperturbed regDB 
#====================================
act.nopert <- as.matrix(read.csv(file='./act/act.regdb.0.csv', row.names = 1)) 
dim(act.nopert)

# Load racipe activities 
act.racipe.df <- read.table(file = '../../../sim.racipe/sim.racipe.comb.kd9/net30tf.states.txt', header = TRUE) 
#rownames(act.racipe.df) <- paste('M', act.racipe.df$MODEL_NO, sep = '') 
rownames(act.racipe.df) <- act.racipe.df$MODEL_NO
act.racipe <- t(act.racipe.df[, rownames(act.nopert)]) # retain TF activities only

cor.nopert <- diag(cor(t(act.racipe), t(act.nopert), method = 'spearman'))
cor.avg.nopert <- mean(abs(diag(cor(t(act.racipe), t(act.nopert), method = 'spearman'))))

saveRDS(cor.nopert, file = paste(outdir, 'cor.nopert.rds', sep=''))
write.csv(format(cor.avg.nopert, digits = NO.SIG.DIGITS), file = paste(outdir, 'avg.cor.nopert.csv', sep=''), row.names = FALSE, quote = FALSE) 

# Calculate cor for all perturbed regulon DBs
#============================================
act.list.by.regdb <- readRDS(file = './act/act.all.rds') # Load activities for all regulon DBs

# First, calculate cor between racipe acvities and calculated activities for each TF and 
# then calculate mean of their abs values
cor.avg.list <- lapply(names(act.list.by.regdb), function(regdb.name) { # go over five regDB names of perturbation levels
  sapply(names(act.list.by.regdb[[regdb.name]]), function (regdb.id) # go over 100 regDBs at each perturbation level 
    mean(abs(diag(cor(t(act.racipe), t(act.list.by.regdb[[regdb.name]][[regdb.id]]), method = 'spearman'))))
  )
}
)
names(cor.avg.list) <- names(act.list.by.regdb)

saveRDS(cor.avg.list, file = paste(outdir, 'avg.cor.pert.rds', sep = ''))

# Calculate median correlations for correlations at each perturbation level
#==========================================================================
medianCor <- sapply(names(cor.avg.list), function(regdb.name) median(cor.avg.list[[regdb.name]])) 
names(medianCor) <- names(cor.avg.list)
medianCor <- cbind(names(medianCor), as.numeric(medianCor))
colnames(medianCor) <- c('pert_level', 'median_cor') 
write.csv(format(medianCor, digits = NO.SIG.DIGITS), file = paste(outdir, 'median.cor.pert.csv', sep=''), row.names = FALSE, quote = FALSE) 
