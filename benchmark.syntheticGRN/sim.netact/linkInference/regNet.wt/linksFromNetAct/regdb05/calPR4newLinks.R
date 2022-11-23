rm(list = ls()) 

NO_BREAKS <- 10 

figdir <- "./figs/"
dir.create(figdir)
outdir <- './data/'

source('./functions.R')

regdb.no <- readRDS(file = '../../../../../regDBs/regDB.nopert.rds')
regdb.all <- readRDS(file = '../../../../../regDBs/regDBs.5.rds')
regdb.per <- regdb.all$regdb.1

# import regulatory network obtained from perturbed regulon
regNet.pert <- construct_regnet_from_regdb(regdb=regdb.per)
regNet.pert <- format_source_target_nodes(inputNet.df = regNet.pert) 
dim(regNet.pert)

# Construct regNet from unperturbed regulons
regNet.df <- construct_regnet_from_regdb(regdb=regdb.no)
dim(regNet.df)
regNet.df <- format_source_target_nodes(inputNet.df = regNet.df) 
dupStatus <- duplicated(regNet.df) 
sum(dupStatus)
regNet.df <- regNet.df[!dupStatus, ] # remove duplicated interactions
dim(regNet.df)
nodes.all.regNet <- sort(unique(union(regNet.df$SOURCE, regNet.df$TARGET)))
length(nodes.all.regNet)  # 447 - nodes are TFs and targets in the regulatory network

# Extract the fraction of regdb (unperturbed) that was replaced for constructing perturbed regdb
regdbList.no <- regdb2String (regdb.cur=regdb.no) # create a list of tf and target pairs
regdbList.per <- regdb2String (regdb.cur=regdb.per) # create a list of tf and target pairs

regdbSet.no <- as.character(unlist(regdbList.no)) # convert the list of tf and target pairs into a set of strings
regdbSet.per <- as.character(unlist(regdbList.per))

# separate the tf and target pairs that are in unperturbed regulon but not in perturbed regulon:
pairsNOPonly <- setdiff(regdbSet.no, regdbSet.per) 
regdbNOPonly <- pairs2regdbFormat(tfTargetPairs=pairsNOPonly)

if(identical(pairsNOPonly, character(0))) {
  print('two regulons are the same')
  quit()
}

# # part of WT regNet Not Perturbed (NOP):
regNetNOPonly <- construct_regnet_from_regdb(regdb=regdbNOPonly) 

dim(regNetNOPonly)


# Calculate precision recall   
#--------------------------- 
inferredRel.df <- read.table(file = paste('./data/inferredLinks.tsv', sep = ''), header = TRUE)
dim(inferredRel.df)


# Keep only thgose links in the inferred links that are NOt in the perturbed regulons 
new_predictions.df <- create_new_predictions (inferredRel.df, regNet.pert)


#inferredRel.df <- inferredRel.df[inferredRel.df$mi>0, ]

# Calculate precision recall  
#--------------------------- 
# Find MI Threshold values
h <- hist(new_predictions.df$mi, breaks = NO_BREAKS, plot=F)  
miTSH.list <- h$breaks[1:(length(h$breaks)-1)]
stats.df <- calPrecisionRecall (inferredRel.df=new_predictions.df, regNet.df=regNetNOPonly, miTSH.list=miTSH.list)
write.table(format(stats.df, digits = 4), file = paste(outdir, 'stats4newLinks.tsv', sep = ''), 
            row.names=TRUE, sep="\t", quote = F, col.names = T)


WIDTH <- 6
HEIGHT <- 6
figname <- paste(figdir, 'stats4newLinks-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')

pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(2,2))
par(mar=c(2.5, 3.0, 1.5, 1.0))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space

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

# precision0-recall curve
plot(stats.df$rCall, stats.df$ppv,
     lty=1, pch=19,
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = 'recall', ylab = 'precision')
lines(stats.df$rCall, stats.df$ppv, lty=2)

hist(new_predictions.df$mi, breaks = NO_BREAKS)  

dev.off()

