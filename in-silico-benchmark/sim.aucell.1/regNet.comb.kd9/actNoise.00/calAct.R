
"Calculate activities using unperturbed and perturbed regulon DBs" 

rm(list = ls()) 
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


outdir <- "./act/"
dir.create(outdir)
actByRegDB_list  <- list()

# RACIPE expressions
#------------------- 
exp_df <- read.table(file = '../../../sim.racipe/sim.racipe.comb.kd9/net30tf.states.exp.txt' , header=TRUE)
obs_ids <- c(paste0('M', exp_df$MODEL_NO)); row.names(exp_df) <- obs_ids # assign row names
exp_df$NO_STATES <- NULL; exp_df$STATE_NO <- NULL ; exp_df$MODEL_NO <- NULL # delete meta columns
exp_df <- t(exp_df) # put genes in rows, models in columns

# Using regulon DBs with NO perturbed targets
#------------------------------------------- 
regdb <- readRDS(file = '../../../regDBsAUCell/regDB.nopert.rds')
act <- calAUcellact(exp_df, regdb)
write.csv(act, file=paste0(outdir, 'act.regdb.0.csv') , quote = FALSE)

# Using regulon DBs with 5 perturbed targets
#------------------------------------------- 
regdb.list <- readRDS(file = '../../../regDBsAUCell/regDBs.5.rds')
actRegDB <- lapply(names(regdb.list), function(regdb_name) calAUcellact(exp_df, regdb.list[[regdb_name]])) 
names(actRegDB) <- names(regdb.list)
saveRDS(actRegDB, file = paste(outdir, 'act.regdb.5.rds', sep = ''))
actByRegDB_list[['regdb.5']] <- actRegDB


# Using regulon DBs with 10 perturbed targets
#------------------------------------------- 
regdb.list <- readRDS(file = '../../../regDBsAUCell/regDBs.10.rds')
actRegDB <- lapply(names(regdb.list), function(regdb_name) calAUcellact(exp_df, regdb.list[[regdb_name]])) 
names(actRegDB) <- names(regdb.list)
saveRDS(actRegDB, file = paste(outdir, 'act.regdb.10.rds', sep = ''))
actByRegDB_list[['regdb.10']] <- actRegDB


# Using regulon DBs with 15 perturbed targets
#------------------------------------------- 
regdb.list <- readRDS(file = '../../../regDBsAUCell/regDBs.15.rds')
actRegDB <- lapply(names(regdb.list), function(regdb_name) calAUcellact(exp_df, regdb.list[[regdb_name]])) 
names(actRegDB) <- names(regdb.list)
saveRDS(actRegDB, file = paste(outdir, 'act.regdb.15.rds', sep = ''))
actByRegDB_list[['regdb.15']] <- actRegDB

saveRDS(actByRegDB_list, file = paste(outdir, 'act.all.rds', sep = ''))
