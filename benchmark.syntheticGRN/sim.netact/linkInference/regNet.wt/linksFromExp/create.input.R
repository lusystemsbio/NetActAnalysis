rm(list = ls()) 

outdir <- './data/'
dir.create(outdir)

# Load and format racipe expressions
exp.racipe <- read.table(file = '../../../../sim.racipe/sim.racipe.wt/net30tf.states.exp.txt', header=TRUE)
mydata <- t(exp.racipe[,4:ncol(exp.racipe)]) # remove metadata columns

write.table(mydata, file = paste(outdir, 'net30tf.states.tsv', sep = ''), row.names=T, sep="\t", quote = F)

