
"Calculate activities using unperturbed and perturbed regulon DBs" 

rm(list = ls()) 
library(viper)
library(dorothea)
source('../../sim.viper/src/regDB2df.R')

outdir <- "./act/"
dir.create(outdir)

# RACIPE expressions
#------------------- 
exp_df <- read.table(file = '../../../sim.racipe/sim.racipe.comb.kd9/net30tf.states.exp.txt' , header=TRUE)
obs_ids <- c(paste0('M', exp_df$MODEL_NO)); row.names(exp_df) <- obs_ids # assign row names
exp_df$NO_STATES <- NULL; exp_df$STATE_NO <- NULL ; exp_df$MODEL_NO <- NULL # delete meta columns
exp_df <- t(exp_df) # put genes in rows, models in columns

# Using regulon DBs with NO perturbed targets
#------------------------------------------- 
regdb <- readRDS(file = '../../../regDBs/regDB.nopert.rds')

# Calculate activities
regulons.df <- regDB2df(regdb = regdb) 
regulons.viper <- df2regulon(regulons.df) 
# act <- aREA(exp_df, regulons.viper) 
# act.es <- act$es
act <- viper(exp_df, regulons.viper, method = "scale", minsize = 20, 
              eset.filter = FALSE, cores = 1, verbose = FALSE)
fname_out <- paste0(outdir, 'act.regdb.0.csv') 
write.csv(act, file=fname_out, quote = FALSE) 

#quit()

# Using regulon DBs with 5 perturbed targets
#------------------------------------------- 
regdb.list <- readRDS(file = '../../../regDBs/regDBs.5.rds')

actRegDB5 <- lapply(names(regdb.list), function(regdb_name) { 
  regulons.df <- regDB2df(regdb = regdb.list[[regdb_name]]); regulons.viper <- df2regulon(regulons.df); 
  #act <- aREA(exp_df, regulons.viper); act.es <- t(act$es)
  act <- viper(exp_df, regulons.viper, method = "scale", minsize = 20, 
                eset.filter = FALSE, cores = 1, 
                verbose = FALSE)
  return(act)
  }
)
names(actRegDB5) <- names(regdb.list)


# Using regulon DBs with 10 perturbed targets
#------------------------------------------- 
regdb.list <- readRDS(file = '../../../regDBs/regDBs.10.rds')

actRegDB10 <- lapply(names(regdb.list), function(regdb_name) { 
  regulons.df <- regDB2df(regdb = regdb.list[[regdb_name]]); regulons.viper <- df2regulon(regulons.df); 
  #act <- aREA(exp_df, regulons.viper); act.es <- t(act$es)
  act <- viper(exp_df, regulons.viper, method = "scale", minsize = 20, 
                eset.filter = FALSE, cores = 1, 
                verbose = FALSE)
  return(act)
}
)
names(actRegDB10) <- names(regdb.list)


# Using regulon DBs with 15 perturbed targets
#------------------------------------------- 
regdb.list <- readRDS(file = '../../../regDBs/regDBs.15.rds')

actRegDB15 <- lapply(names(regdb.list), function(regdb_name) { 
  regulons.df <- regDB2df(regdb = regdb.list[[regdb_name]]); regulons.viper <- df2regulon(regulons.df); 
  #act <- aREA(exp_df, regulons.viper); act.es <- t(act$es)
  act <- viper(exp_df, regulons.viper, method = "scale", minsize = 20, 
                eset.filter = FALSE, cores = 1, 
                verbose = FALSE)
  return(act)
}
)
names(actRegDB15) <- names(regdb.list)

# Save activities
#----------------
saveRDS(actRegDB5, file = paste(outdir, 'act.regdb.5.rds', sep = ''))
saveRDS(actRegDB10, file = paste(outdir, 'act.regdb.10.rds', sep = ''))
saveRDS(actRegDB15, file = paste(outdir, 'act.regdb.15.rds', sep = '')) 

actByRegDB_list  <- list()
actByRegDB_list[['regdb.5']] <- actRegDB5
actByRegDB_list[['regdb.10']] <- actRegDB10
actByRegDB_list[['regdb.15']] <- actRegDB15
saveRDS(actByRegDB_list, file = paste(outdir, 'act.all.rds', sep = ''))
