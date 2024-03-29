rm(list = ls()) 

NO_BREAKS <- 10

figdir <- "./figs/"
dir.create(figdir)
outdir <- './data/'

source('./functions.R')

# Load unperturbed regulons and construct regulatory links
regdb.no <- readRDS(file = '../../../../../regDBs/regDB.nopert.rds')
regNet.df <- construct_regnet_from_regdb(regdb=regdb.no)
dim(regNet.df)
regNet.df <- format_source_target_nodes(inputNet.df = regNet.df) 

dupStatus <- duplicated(regNet.df)
sum(dupStatus) 
regNet.df <- regNet.df[!dupStatus, ] # remove duplicated interactions
dim(regNet.df)
nodes.all.regNet <- sort(unique(union(regNet.df$SOURCE, regNet.df$TARGET)))
length(nodes.all.regNet)  # 447 - nodes are TFs and targets in the regulatory network

# Calculate precision recall   
#--------------------------- 
inferredRel.df <- read.table(file = paste('./data/inferredLinks.tsv', sep = ''), header = TRUE)
dim(inferredRel.df)

#inferredRel.df <- inferredRel.df[inferredRel.df$mi>0, ]

# Calculate precision recall  
#--------------------------- 
# Find MI Threshold values
h <- hist(inferredRel.df$mi, breaks = NO_BREAKS, plot=FALSE)  
miTSH.list <- h$breaks[1:(length(h$breaks)-1)]

ppvRcall <- calPrecisionRecall(inferredRel.df=inferredRel.df, regNet.df=regNet.df, miTSH.list=miTSH.list)

write.table(format(ppvRcall, digits = 4), file = paste(outdir, 'prValues.tsv', sep = ''), 
            row.names=TRUE, sep="\t", quote = F, col.names = T)

WIDTH <- 6
HEIGHT <- 6
figname <- paste(figdir, 'stats-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')

pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(2,2))
par(mar=c(2.5, 3.0, 1.5, 1.0))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space

# precision0-recall curve
plot(ppvRcall$rCall, ppvRcall$ppv,
     lty=1, pch=19,
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = 'recall', ylab = 'precision')
lines(ppvRcall$rCall, ppvRcall$ppv, lty=2)

hist(inferredRel.df$mi, main = '')

performVal <- ppvRcall$rCall
plot(performVal,
     lty=1, pch=19,
     ylim = c(0, 1), xlab = 'mutual information threshold', ylab = 'recall',
     xaxt = "n")
axis(1,                         # Define x-axis manually
     at = 1:length(miTSH.list),
     labels = miTSH.list)
lines(performVal, lty=2)

performVal <- ppvRcall$ppv
plot(performVal, #xlim = c(1, length(rCall.list)),
     lty=1, pch=19,
     ylim = c(0, 1), xlab = 'mutual information threshold', ylab = 'precision',
     xaxt = "n")
axis(1,                         # Define x-axis manually
     at = 1:length(miTSH.list),
     labels = miTSH.list)
lines(performVal, lty=2)

dev.off()

