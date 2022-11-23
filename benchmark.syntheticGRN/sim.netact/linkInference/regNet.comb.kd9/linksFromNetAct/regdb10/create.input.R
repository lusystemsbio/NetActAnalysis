rm(list = ls()) 
outdir <- './data/'
dir.create(outdir)

# NetAct calculated activities are used as TF activities and 
# target expressions are used as their activities 

# Load racipe expressions
exp.racipe <- read.table(file = '../../../../../sim.racipe/sim.racipe.comb.kd9/net30tf.states.exp.txt', header=TRUE)
exp.racipe <- exp.racipe[,4:ncol(exp.racipe)] # remove metadata columns
dim(exp.racipe)

# Load TF activities at specific perturbation level
act.regdb.all <- readRDS(file = '../../../../regNet.comb.kd9/actNoise.00/act/act.regdb.10.rds') 
act.tf <- act.regdb.all[["regdb.1"]]
act.tf <- t(act.tf)

# replace TF expression with TF activity 
exp.targets <- exp.racipe[,setdiff(colnames(exp.racipe), colnames(act.tf))]
dim(exp.targets)

# Combine  TF activity and target expressions
mydata <- t(cbind(act.tf, exp.targets))

write.table(mydata, file = paste(outdir, 'net30tf.states.tsv', sep = ''), row.names=T, sep="\t", quote = F)

# Extract TF list 
#regdb <- readRDS('../../../../regDBs/regDB.nopert.rds')
#names(regdb)

# Save only TFs
#tfs <- names(regdb)
#write.table(tfs, file = paste(outdir, 'transcription_factors.tsv', sep = ''), 
#            row.names=FALSE, sep="\t", quote = F, col.names = F)

