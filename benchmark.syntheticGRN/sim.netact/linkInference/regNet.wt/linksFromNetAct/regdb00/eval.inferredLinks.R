rm(list = ls()) 

NO_BREAKS <- 10

figdir <- "./figs/"
dir.create(figdir)
outdir <- './data/'

source('./functions.R')

# Load unperturbed regulons and construct regulatory links
regdb.pos <- readRDS(file = '../../../../../regDBs/regDB.nopert.rds')
regdbList.pos <- regdb2String (regdb.cur=regdb.pos) # create a list of tf and target pairs
regdbSet.pos <- as.character(unlist(regdbList.pos)) # convert the list of tf and target pairs into a set of strings

# find negative regulons
targets.all <- unique(union(names(regdb.pos), unlist(regdb.pos)))
regdb.neg <- sapply(names(regdb.pos), function(tf) setdiff(targets.all, regdb.pos[[tf]]))
regdbList.neg <- regdb2String (regdb.cur=regdb.neg)
regdbSet.neg <- as.character(unlist(regdbList.neg))
length(regdbSet.neg)

gtSet <- list()
gtSet[['regdbSet.pos']] <- regdbSet.pos
gtSet[['regdbSet.neg']] <- regdbSet.neg

# Load inferred links
#-------------------- 
inferredRel.df <- read.table(file = paste('./data/inferredLinks.tsv', sep = ''), header = TRUE)
dim(inferredRel.df)

#inferredRel.df <- inferredRel.df[inferredRel.df$mi>0, ] 

# Calculate precision recall  
#--------------------------- 
# Find MI Threshold values
h <- hist(inferredRel.df$mi, breaks = NO_BREAKS, plot=F)  
miTSH.list <- h$breaks[1:(length(h$breaks)-1)]


#stats.df <- evaluatePredictions(inferredRel.df=inferredRel.df, gtLinks.df=regNet.df, miTSH.list=miTSH.list)
stats.df <- evaluatePredictions(inferredRel.df=inferredRel.df, gtSet=gtSet, miTSH.list=miTSH.list)

write.table(format(stats.df, digits = 4), file = paste(outdir, 'stats.tsv', sep = ''), 
            row.names=TRUE, sep="\t", quote = F, col.names = T)



WIDTH <- 4
HEIGHT <- 8
figname <- paste(figdir, 'stats-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')

pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(2,1))
par(mar=c(2.5, 3.0, 1.5, 1.0))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space

# precision0-recall curve
plot(stats.df$rCall, stats.df$ppv,
     lty=1, pch=19,
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = 'recall', ylab = 'precision')
lines(stats.df$rCall, stats.df$ppv, lty=2)

# ROC curve
plot(stats.df$fpr, stats.df$tpr,
     lty=1, pch=19,
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = 'False positive rate', ylab = 'True positive rate')
lines(stats.df$fpr, stats.df$tpr, lty=2)
dev.off()




WIDTH <- 6
HEIGHT <- 6
figname <- paste(figdir, 'prCurve.gt-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')

pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(2,2))
par(mar=c(2.5, 3.0, 1.5, 1.0))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space

# precision0-recall curve
plot(stats.df$rCall, stats.df$ppv,
     lty=1, pch=19,
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = 'recall', ylab = 'precision')
lines(stats.df$rCall, stats.df$ppv, lty=2)

hist(inferredRel.df$mi, main = '')

performVal <- stats.df$rCall
plot(performVal,
     lty=1, pch=19,
     ylim = c(0, 1), xlab = 'mutual information threshold', ylab = 'recall',
     xaxt = "n")
axis(1,                         # Define x-axis manually
     at = 1:length(miTSH.list),
     labels = miTSH.list)
lines(performVal, lty=2)

performVal <- stats.df$ppv
plot(performVal, #xlim = c(1, length(rCall.list)),
     lty=1, pch=19,
     ylim = c(0, 1), xlab = 'mutual information threshold', ylab = 'precision',
     xaxt = "n")
axis(1,                         # Define x-axis manually
     at = 1:length(miTSH.list),
     labels = miTSH.list)
lines(performVal, lty=2)

dev.off()
