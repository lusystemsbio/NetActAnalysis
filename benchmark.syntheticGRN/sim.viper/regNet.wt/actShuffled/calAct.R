
"Calculate activities using unperturbed regulon DB for racipe expression by shuffling gene labels" 

rm(list = ls()) 
library(viper) 
library(dorothea) 
source('../../src/regDB2df.R')


NUM_PERMUTATIONS <- 10000
#NUM_PERMUTATIONS <- 100 

outdir <- "./act/" 
dir.create(outdir) 

set.seed(100)

# RACIPE expressions
#------------------- 
exp_df <- read.table(file = '../../../sim.racipe/sim.racipe.wt/net30tf.states.exp.txt', header=TRUE)
obs_ids <- c(paste0('M', exp_df$MODEL_NO)); row.names(exp_df) <- obs_ids # assign row names
exp_df$NO_STATES <- NULL; exp_df$STATE_NO <- NULL; exp_df$MODEL_NO <- NULL # delete meta columns
exp_df <- t(exp_df) # put genes in rows, models in columns

# Using regulon DBs with NO perturbed targets
#------------------------------------------- 
regdb <- readRDS(file = '../../../regDBs/regDB.nopert.rds') # Load Regulon DB

regulons.df <- regDB2df(regdb = regdb) 
regulons.viper <- df2regulon(regulons.df) 

# calculate TF activities 
gene_names <- rownames(exp_df)
act_list <- lapply(1:NUM_PERMUTATIONS, function(count_perm) {
  rownames(exp_df) <- sample(gene_names)
  act <- viper(exp_df, regulons.viper, method = "scale", minsize = 20, 
                eset.filter = FALSE, cores = 1, verbose = FALSE)
  return(act)
}
)
names(act_list) <- 1:NUM_PERMUTATIONS

outdir <- "./act/" 
dir.create(outdir)
saveRDS(act_list, file = paste(outdir, 'act.all.rds', sep = '')) 
