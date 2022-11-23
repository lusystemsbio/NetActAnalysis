"
   Calculate activities for unperturbed and perturbed regulon DBs
"
rm(list=ls()) 

NO_PERMUTATIONS <- 1000
set.seed(100)

source('./lib.tfnetwork.dist.R')

outdir <- './act/' 
dir.create(outdir) 

act.all  <- list() # List to save activities for each reg DB

# Load reg DB 0 
regdb <- readRDS(file = '../../../regDBs/regDB.nopert.rds')
length(names(regdb))
length(unique(unlist(regdb)))

# Create and save connectivity matrix (one time):
ret_list <- create_conn_matrix(regdb)
conn_df <- ret_list[["CONN.MAT"]]
targetGenes <- ret_list[["TARGET.LIST"]]
write.csv(conn_df, file='./net30tf.matrix.csv', quote = FALSE)

# Load RACIPE EXPR
exp_df <- read.table(file = '../../../sim.racipe/sim.racipe.wt/net30tf.states.exp.txt', header=TRUE)
row.names(exp_df) <- c(paste0('M', exp_df$MODEL_NO)) 
exp_df$NO_STATES <- NULL; exp_df$STATE_NO <- NULL ; exp_df$MODEL_NO <- NULL # delete meta columns
exp_df <- t(exp_df[,targetGenes]) # put genes in rows, models in columns

# Calculate activities
act.all <- vector(mode = 'list', length = NO_PERMUTATIONS) 

for (k in seq(1, NO_PERMUTATIONS)) {
   # shuffle gene labels:
   targetGenesShuffled <- sample(targetGenes)
   exp_df <- exp_df[targetGenesShuffled, ] 
   rownames(exp_df) <- targetGenesShuffled  
   
   # save the expressions in the NCA format (columns are models)
   write.csv(exp_df, file='./net30tf.states.exp.csv', quote = FALSE) 
   
   # calculate NCA activities:
   system("octave calFastNCA.m") 
   
   act.all[[k]] <- format.NCAoutput()  
}
saveRDS(act.all, paste(outdir, 'act.all.rds', sep = '')) 
