"Calculate Activites for the TFs between racipe and netact expressions"

rm(list = ls())
library(NetAct)

NUM_PERMUTATIONS <- 1000

set.seed(100) 

outdir <- "./act/"
dir.create(outdir)

# Load RACIPE expressions
#------------------------- 
exp_df <- read.table(file = '../../../sim.racipe/sim.racipe.wt/net30tf.states.exp.txt', header=TRUE)
obs_ids <- c(paste0('M', exp_df$MODEL_NO)); row.names(exp_df) <- obs_ids # create row names
exp_df$NO_STATES <- NULL; exp_df$STATE_NO <- NULL ; exp_df$MODEL_NO <- NULL # delete meta columns
exp_df <- t(exp_df) # put genes in rows, models in columns

# create DEresult object 
#-----------------------
# Create an DE result list object containing a table with on column named 'padj'
DEresult <- list()
myTable <- as.data.frame(matrix(nrow = nrow(exp_df), ncol = 1), 
                         col.names=c('padj'))
colnames(myTable) <- c('padj') 
rownames(myTable) <- rownames(exp_df)
myTable$padj <- 1 
DEresult$table <- myTable

# Load regulon DB
regdb <- readRDS(file = '../../../regDBs/regDB.nopert.rds')

# Calculate TF activity using gene permutated expression data
gene_names <- rownames(exp_df)
act.all <- vector(mode = "list", length = NUM_PERMUTATIONS)
act.all <- lapply(1:NUM_PERMUTATIONS, function(count_perm) {
  rownames(exp_df) <- sample(gene_names)
  acts = TF_Activity(names(regdb), regdb, exp_df, DEresult, with_weight = FALSE)
  return(acts$all_activities)
 }
)
names(act.all) <- 1:NUM_PERMUTATIONS

saveRDS(act.all, file = paste(outdir, 'act.all.rds', sep = ''))
