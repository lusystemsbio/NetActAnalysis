"Calculate Activites for the TFs between racipe and netact expressions"

rm(list = ls())
# library(NetAct)

NUM_PERMUTATIONS <- 1000

regDBdir <- '../../../regDBsAUCell/'


library(AUCell)

calAUcellact <- function(exp_df, regdb){
  cells_rankings <- AUCell_buildRankings(exp_df, plotStats = F)
  class(cells_rankings)
  length(cells_rankings)
  rankings <- getRanking(cells_rankings) 
  
  aucMaxRank <- ceiling(nrow(rankings)*5/100) # 5%
  aucMaxRank
  cells_AUC <- AUCell_calcAUC(regdb, cells_rankings, aucMaxRank=aucMaxRank, nCores=1) 
  act <- getAUC(cells_AUC)
  return(act)
}


set.seed(100) 

outdir <- "./act/"
dir.create(outdir)

# Load RACIPE expressions
#------------------------- 
exp_df <- read.table(file = '../../../sim.racipe/sim.racipe.wt/net30tf.states.exp.txt' , header=TRUE)
obs_ids <- c(paste0('M', exp_df$MODEL_NO)); row.names(exp_df) <- obs_ids # create row names
exp_df$NO_STATES <- NULL; exp_df$STATE_NO <- NULL ; exp_df$MODEL_NO <- NULL # delete meta columns
exp_df <- t(exp_df) # put genes in rows, models in columns

# Load regulon DB
regdb <- readRDS(file = paste0(regDBdir, 'regDB.nopert.rds'))

# Calculate TF activity using gene permutated expression data
gene_names <- rownames(exp_df)
act.all <- vector(mode = "list", length = NUM_PERMUTATIONS)
act.all <- lapply(1:NUM_PERMUTATIONS, function(count_perm) {
  rownames(exp_df) <- sample(gene_names)
  calAUcellact(exp_df, regdb)
 }
)
names(act.all) <- 1:NUM_PERMUTATIONS
saveRDS(act.all, file = paste(outdir, 'act.all.rds', sep = '')) 
