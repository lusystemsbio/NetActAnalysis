
rm(list = ls())

outdir <- "./data/" 
dir.create(outdir) 

# Load racipe expressions
mydata <- read.table(file = '../../../sim.racipe/sim.racipe.wt/net30tf.states.exp.txt', header=TRUE)
mydata<-mydata[,4:ncol(mydata)]

write.table(mydata, file = paste(outdir, 'net30tf.states.tsv', sep = ''), row.names=FALSE, sep="\t", quote = F)

# load regulon DB 
regdb <- readRDS('../../../regDBs/regDB.nopert.rds')
names(regdb)

# Save only TFs
write.table(names(regdb), file = paste(outdir, 'transcription_factors.tsv', sep = ''), 
             row.names=FALSE, sep="\t", quote = F, col.names = F)

