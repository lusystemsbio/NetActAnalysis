
rm(list = ls())

NO_BREAKS <- 6

source('./functions.R')

figdir <- './figs/'
dir.create(figdir)
outdir <- './data/'

regNetDir <- '../../../../regNet/data/'
dir.create(outdir) 

inferredNet.df <- read.delim(file = './output.tsv', sep = '\t')  
inferredNet.df <- inferredNet.df[, c(2, 3, 4)] # first column 
inferredNet.df <- format_source_target_nodes(inputNet.df = inferredNet.df)  
inferredNet.df$importance <- as.numeric(inferredNet.df$importance)
dupStatus <- duplicated(inferredNet.df )
sum(dupStatus)
dim(inferredNet.df)
inferredNet.df <- inferredNet.df[!dupStatus, ] # remove duplicated interactions
dim(inferredNet.df)
nodes.all <- unique(union(inferredNet.df$TF, inferredNet.df$target))
nodes.all <- sort(nodes.all)
length(nodes.all) # 30 - only TFs


regNet.df <- read.delim2(file = paste0(regNetDir, 'net30tf.tpo'), sep = '\t')
regNet.df <- format_source_target_nodes(inputNet.df = regNet.df) 
regNet.df <- regNet.df[,c(1,2)] # retain only source and target
dupStatus <- duplicated(regNet.df)
sum(dupStatus)
regNet.df <- regNet.df[!dupStatus, ] # remove duplicated interactions
dim(regNet.df)
nodes.all.regNet <- sort(unique(union(regNet.df$SOURCE, regNet.df$TARGET)))
length(nodes.all.regNet)  # 447 - nodes are TFs and targets in the regulatory network

# remove interactions involving TF targets (starting with tg)
# regNet.df <- regNet.df[!(startsWith(regNet.df$SOURCE, 'tg') | startsWith(regNet.df$TARGET, 'tg')), ]
# dim(regNet.df)

inferredNet.df$importance <- (inferredNet.df$importance-min(inferredNet.df$importance))/(max(inferredNet.df$importance)-min(inferredNet.df$importance))

h <- hist(inferredNet.df$importance, breaks = NO_BREAKS, plot=F) 
length(h$breaks)

TSH.list <- h$breaks[1:(length(h$breaks)-1)]

# calculate precision and recall
ppvRcall <- sapply(TSH.list, function(tsh) {
  #tsh <- TSH.list[1]
  inferredNet.retained.df <- inferredNet.df[inferredNet.df$importance>tsh, c(1, 2)] 
  colnames(inferredNet.retained.df) <- colnames(regNet.df) 
  dupStatus <- duplicated(inferredNet.retained.df)
  inferredNet.retained.df <- inferredNet.retained.df[!dupStatus,] # remove duplicated interactions
  
  combInt.df <- rbind(regNet.df, inferredNet.retained.df)
  dupStatus <- duplicated(combInt.df)
  if(sum(dupStatus) > 0){
    rCall <- sum(dupStatus)/nrow(regNet.df) # calculate recall 
    ppv <- sum(dupStatus)/nrow(inferredNet.retained.df) # calculate precision/positive predictive value
    ppvRcall <- c("ppv"=ppv, "rCall"=rCall)
  } else {
    rCall <- ppv <- NA 
  }
  ppvRcall <- c("ppv"=ppv, "rCall"=rCall)
  return(ppvRcall)
})  

class(ppvRcall)
dim(ppvRcall)
ppvRcall <- t(ppvRcall)
class(ppvRcall)
class(ppvRcall[,'rCall'])
rownames(ppvRcall) <- TSH.list


write.table(format(ppvRcall, digits = 4), file = paste(outdir, 'precision-recall.tsv', sep = ''), 
            row.names=TRUE, sep="\t", quote = F, col.names = T)



WIDTH <- 8
HEIGHT <- 8
figname <- paste(figdir, 'stats-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width=WIDTH , height=HEIGHT, paper='special', onefile = TRUE)
par(mfrow=c(2,2)) 
par(mar=c(2.5, 3.0, 1.5, 1.0))  # bottom, left, top, right
par(mgp=c(1.5, 0.4, 0))
par(oma=c(3, 1, 1, 1)) # all sides have 3 lines of space

# Precision-Recall curve
plot(ppvRcall[,'rCall'], ppvRcall[,'ppv'], 
     lty=1, pch=19,
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = 'recall', ylab = 'precision')
lines(ppvRcall[,'rCall'], ppvRcall[,'ppv'], lty=2) 

#hist(inferredNet.df$importance, xlab = 'importance score', ylab = 'count' , main = '')
hist(inferredNet.df$importance, xlab = 'normed importance score', ylab = 'count' , main = '')

# Precision 
plot(ppvRcall[,'ppv'], 
     lty=1, pch=19, 
     ylim = c(0, 1), xlab = 'importance score threshold', ylab = 'Precision', 
     xaxt = "n")  
axis(1,                         # Define x-axis manually
     at = 1:length(TSH.list),
     labels = TSH.list)
lines(ppvRcall[,'ppv'], lty=2) 

plot(ppvRcall[,'rCall'],  
     lty=1, pch=19, 
     ylim = c(0, 1), xlab = 'importance score threshold', ylab = 'recall', 
     xaxt = "n")  
axis(1,                         # Define x-axis manually
     at = 1:length(TSH.list),
     labels = TSH.list)
lines(ppvRcall[,'rCall'], lty=2)

dev.off() 

