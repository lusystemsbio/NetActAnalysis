rm(list = ls()) 
outdir <- './data/'
dir.create(outdir)

# NetAct calculated activities are used as TF activities and 
# target expressions are used as their activities 
# Load racipe expressions
exp.racipe <- read.table(file = '../../../../../sim.racipe/sim.racipe.wt/net30tf.states.exp.txt', 
                         header=TRUE)
exp.racipe <- exp.racipe[,4:ncol(exp.racipe)] # remove metadata columns
dim(exp.racipe)

# Load TF  activities  
act.tf <- read.csv(file = '../../../../regNet.wt/actNoise.00/act/act.regdb.0.csv', 
                   header=TRUE, row.names = 1)
class(act.tf )
act.tf  <- t(act.tf)
dim(act.tf)

# replace TF expression with TF activity 
exp.targets <- exp.racipe[,setdiff(colnames(exp.racipe), colnames(act.tf))]
dim(exp.targets)

#sort(colnames(exp.targets))
head(sort(colnames(exp.targets)))
tail(sort(colnames(exp.targets)))

# Combine  TF activity and target expressions
mydata <- t(cbind(act.tf, exp.targets))

write.table(mydata, file = paste(outdir, 'net30tf.states.tsv', sep = ''), row.names=T, sep="\t", quote = F)

# Extract TFs
#regdb <- readRDS('../../../../../regDBs/regDB.nopert.rds')
#tfs <- names(regdb)
#write.table(tfs, file = paste(outdir, 'transcription_factors.tsv', sep = ''), 
#            row.names=FALSE, sep="\t", quote = F, col.names = F)
