rm(list=ls()) 

#set.seed(100)

outdir <- './data/'
dir.create(outdir)

fname_network_csv <- paste0(outdir, 'net30tf.csv') 
fname_network_tpo <- paste0(outdir, 'net30tf.tpo') 
fname_regulon <- paste0(outdir, 'regdb.net30tf.rda') 
fname_network_comp_csv <- paste0(outdir, 'net30tf.comp.csv')  # complement network:

source('./lib.tfnetwork.R')

# CONSTANTS 
NO_TFS <- 30
NO_LINKS_PER_TF <- 2
TFSET <- sprintf("tf%d", 1:NO_TFS)  
TARGET_POOLSIZE <- 500
NO_TARGETS_PER_TF <- 20
TARGETSET <- sprintf("tg%d", 1:TARGET_POOLSIZE) 



# Create  regulatory network and build the corresponding regulon 
tfnet <- create_tfnet(TFSET, NO_LINKS_PER_TF) # create links among TFs
tf2target_net <- create_tf2target_net(TFSET, TARGETSET, NO_TARGETS_PER_TF) # create links from TF to targets
tfnet = rbind(tfnet, tf2target_net) # combine the interactions within TF set with those between TF set and TARGET set
reg_db <- create_regdb (tfnet) # infer regulon DB from the created network 

write.csv(tfnet, file=fname_network_csv, quote = FALSE, row.names = FALSE)
write.table(tfnet, file=fname_network_tpo, quote = FALSE, row.names = FALSE, sep='\t')
 
save(reg_db, file = fname_regulon)
saveRDS(reg_db, file = paste(outdir, "regdb.nopert.rds", sep = ''))

#====================================================================#
#  Create COMPLEMENTARY NETWORK (consisting of non-interacting pairs)
tfnet_comp <- data.frame('SOURCE' = character(), 'TARGET' = character()) 
length(TFSET)
TARGETSET_USED <- unique(tf2target_net$TARGET)
length(TARGETSET_USED)

TARGETSET_USED <- unique(tf2target_net$TARGET)
length(TARGETSET_USED)

for (g in TFSET){ 
   # TF subset not used for LINKS among TFSET
   target_set1 <- setdiff(TFSET, reg_db[[g]])
   
   # target subset not used for LINKS from TFSET to TARGETSET:
   target_set2 <- setdiff(TARGETSET_USED, reg_db[[g]]) 
   
   # merge the two subsets: 
   targets <- union(target_set1, target_set2)
   tfsubnet <- data.frame('SOURCE' = rep(g, each=length(targets)), 
                          'TARGET' = targets)
   tfnet_comp <- rbind(tfnet_comp, tfsubnet)
   #break
}

tfnet_comp$TYPE <- 0 
write.csv(tfnet_comp, file=fname_network_comp_csv, quote = FALSE, row.names = FALSE)
