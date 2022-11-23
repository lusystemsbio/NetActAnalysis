
rm(list = ls()) 

outdir <- "./data/" 
dir.create(outdir) 

# Load racipe expressions
exp.racipe <- read.table(file = '../../../../sim.racipe/sim.racipe.comb.kd9/net30tf.states.exp.txt', header=TRUE)
exp.racipe <- exp.racipe[,4:ncol(exp.racipe)] # remove metadata columns
dim(exp.racipe)

# Load TF activities
act.regdb.all <- readRDS(file = '../../../../sim.netact/regNet.comb.kd9/actNoise.00/act/act.regdb.5.rds') 
act.tf <- act.regdb.all[["regdb.1"]]
act.tf <- t(act.tf)

# Separate target expressions from racipe expressions
exp.targets <- exp.racipe[,setdiff(colnames(exp.racipe), colnames(act.tf))]
dim(exp.targets)

# Combine  TF activity and target expressions
mydata <- cbind(act.tf, exp.targets)

write.table(mydata, file = paste(outdir, 'net30tf.states.tsv', sep = ''), row.names=FALSE, sep="\t", quote = F)

# load regulon DB 
regdb <- readRDS('../../../../regDBs/regDB.nopert.rds')
names(regdb)

write.table(names(regdb), file = paste(outdir, 'transcription_factors.tsv', sep = ''), 
            row.names=FALSE, sep="\t", quote = F, col.names = F)
