rm(list = ls())

NO_BREAKS <- 10

source('./functions.R')

outdir <- inputdir <-  './data/'
regNetDir <- '../../../regNet/data/' 

figdir <- "./figs/" 
dir.create(figdir) 

# Extract regulon db with no perturbation
regdb <- readRDS('../../../regDBs/regDB.nopert.rds')
names(regdb)
 
# Load topology file and keep the only TF-TF interactions 
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
#regNet.df <- regNet.df[!(startsWith(regNet.df$SOURCE, 'tg') | startsWith(regNet.df$TARGET, 'tg')), ]
#dim(regNet.df)


# calculate precision recall using inferred links  
#------------------------------------------------ 
inferredRel.df <- read.table(file = paste(inputdir, 'inferredRel.tsv', sep = ''), header = TRUE)

inferredRel.df$mi <- abs(inferredRel.df$corVal) 

# retain sources that are TFs 
names(regdb)
inferredRel.df <- inferredRel.df[inferredRel.df$Gene1 %in% names(regdb), ] 
dim(inferredRel.df)

min(inferredRel.df$corVal)
max(inferredRel.df$corVal)


# Find MI Threshold values
h <- hist(inferredRel.df$mi, breaks = NO_BREAKS, plot=F)  
max(inferredRel.df$mi) 
min(inferredRel.df$mi)
h$breaks
h$counts
length(h$breaks)
miTSH.list <- h$breaks[1:(length(h$breaks)-1)]
miTSH.list


ppvRcall <- calPrecisionRecall(inferredRel.df=inferredRel.df, regNet.df=regNet.df, miTSH.list=miTSH.list)
write.table(format(ppvRcall, digits = 4), file = paste(outdir, 'precision-recall.tsv', sep = ''), 
            row.names=TRUE, sep="\t", quote = F, col.names = T)



# plot PR curves 
WIDTH <- 6
HEIGHT <- 6
figname <- paste(figdir, 'stats-' ,WIDTH, 'x', HEIGHT,'.pdf', sep = '')

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
     ylim = c(0, 1), xlab = 'abs(partial cor) TSH', ylab = 'recall',
     xaxt = "n")
axis(1,                         # Define x-axis manually
     at = 1:length(miTSH.list),
     labels = miTSH.list)   
lines(performVal, lty=2)

performVal <- ppvRcall$ppv
plot(performVal, #xlim = c(1, length(rCall.list)),
     lty=1, pch=19,
     ylim = c(0, 1), xlab = 'abs(partial cor) TSH', ylab = 'precision',
     xaxt = "n")
axis(1,                         # Define x-axis manually
     at = 1:length(miTSH.list),
     labels = miTSH.list) 
lines(performVal, lty=2)
dev.off()

