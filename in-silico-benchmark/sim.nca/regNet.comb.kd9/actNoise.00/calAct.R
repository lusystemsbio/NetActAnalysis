"
   Calculate activities for unperturbed and perturbed regulon DBs
"
rm(list=ls()) 

source('./lib.tfnetwork.dist.R')

outdir <- './act/' 
dir.create(outdir) 

act.all  <- list() # List to save activities for each reg DB

# Load RACIPE EXPR
exp_df <- read.table(file = '../../../sim.racipe/sim.racipe.comb.kd9/net30tf.states.exp.txt', header=TRUE)
row.names(exp_df) <- c(paste0('M', exp_df$MODEL_NO)) 
exp_df$NO_STATES <- NULL; exp_df$STATE_NO <- NULL ; exp_df$MODEL_NO <- NULL # delete meta columns
exp_df <- t(exp_df) # put genes in rows, models in columns

# Calculate act: reg DB 0 
regdb <- readRDS(file = '../../../regDBs/regDB.nopert.rds')
length(names(regdb))
length(unique(unlist(regdb)))

create.inputs.forFastNCA(regdb=regdb, exp_df=exp_df) 
system("octave calFastNCA.m") # Calculate NCA activities 
act <- format.NCAoutput()
write.csv(act, file= paste0(outdir, 'act.regdb.0.csv'), quote = FALSE)


# Calculate act: reg DB 5
regdb.list <- readRDS(file = '../../../regDBs/regDBs.5.rds')
act.regdb <- cal.nca.act(regdb.list, exp_df)  
saveRDS(act.regdb, paste(outdir, 'act.regdb.5.rds', sep = ''))
act.all[['regdb.5']] <- act.regdb 

#quit()

# Calculate act: reg DB 10
regdb.list <- readRDS(file = '../../../regDBs/regDBs.10.rds')
act.regdb <- cal.nca.act(regdb.list, exp_df) 
saveRDS(act.regdb, paste(outdir, 'act.regdb.10.rds', sep = ''))
act.all[['regdb.10']] <- act.regdb

# Calculate act: reg DB 15
regdb.list <- readRDS(file = '../../../regDBs/regDBs.15.rds')
act.regdb <- cal.nca.act(regdb.list, exp_df) 
saveRDS(act.regdb, paste(outdir, 'act.regdb.15.rds', sep = ''))
act.all[['regdb.15']] <- act.regdb

saveRDS(act.all, paste(outdir, 'act.all.rds', sep = '')) 
