#===========================================================#
#                    LIBRARY FOR FUNTION DEFINITIONS  
#                          lib.tfnetwork.R
#                        tfnetwork/lib.r
#===========================================================#

#---------------------------------------------------------------------#
#                 compute correlation distance 
#---------------------------------------------------------------------#
dist_cor <- function(x) as.dist(1-cor(t(x)))

#---------------------------------------------------------------------#
#        normalize the expressions across (columns)
#---------------------------------------------------------------------#

normalize_across_columns <- function(data_p){
   #normalize data_l and return the normalized data
   ndata_l=data_p
   for(i in 1:dim(data_p)[2]) 
   {
      m <- mean(data_p[,i])
      sdt <- sd(data_p[,i]) 
      ndata_l[,i] = (data_p[,i]-m)/sdt
   } 
   return(ndata_l)
}

#---------------------------------------------------------------------#
#           normalize the expressions across rows
#---------------------------------------------------------------------#
normalize_across_rows <- function(data_p){
   #normalize data_l and return the normalized data
   ndata_l=data_p
   for(i in 1:dim(data_p)[1]) 
   {
      #m <- mean(data_p[,i])
      m=mean(as.numeric(data_p[i,]))
      #sdt <- sd(data_p[,i]) 
      sdt=sd(as.numeric(data_p[i,])) 
      
      #ndata_l[,i] = (data_p[,i]-m)/sdt
      ndata_l[i,] = (data_p[i,]-m)/sdt
   } 
   return(ndata_l)
}

#---------------------------------------------------------------------#
#           create TF network
#---------------------------------------------------------------------#
create_tfnet <- function(TFSET, NO_LINKS_PER_TF) {
   TYPE_ACT <- 1
   TYPE_INH <- 2
   TYPE_SIG <- 5
   
   NO_TFS <- length(TFSET)
   TOTAL_LINKS <- NO_TFS * NO_LINKS_PER_TF
   
   PERCENT_ACT_LINKS <- 0.25 # percent of activation links 1
   PERCENT_INH_LINKS <- 0.25 # percent of inhibition links  2
   PERCENT_SIG_LINKS <- 0.50 # percent of phosphorylation signaling links 5
   
   NO_ACT_LINKS <- TOTAL_LINKS * PERCENT_ACT_LINKS
   NO_INH_LINKS <- TOTAL_LINKS * PERCENT_INH_LINKS
   NO_SIG_LINKS <- TOTAL_LINKS * PERCENT_SIG_LINKS
   

   #TFSET <- sprintf("tf%d", 1:NO_TFS)  
   
   # create empy data frame for TF network: 
   tfnet <- data.frame('SOURCE' = character(), 
                       'TARGET' = character())
   
   # establish interactions among TFs: 
   for (g in TFSET){ 
      targets <- sample(TFSET, NO_LINKS_PER_TF)
      tfsubnet <- data.frame('SOURCE' = rep(g, each=NO_LINKS_PER_TF), 
                             'TARGET' = targets)
      tfnet <- rbind(tfnet, tfsubnet)
      #reg_db[[g]] <- targets
   }
   
   # add a column LINK_ID: 
   tfnet <- cbind(LINK_ID = seq(1:TOTAL_LINKS), tfnet)
   
   # add a column TYPE: 
   tfnet <- cbind(tfnet, TYPE = 0)
   
   act_link_ids <- sample(tfnet$LINK_ID, NO_ACT_LINKS, replace=FALSE) 
   
   remaining_link_ids <- setdiff(tfnet$LINK_ID, act_link_ids)
   inh_link_ids <- sample(remaining_link_ids, NO_INH_LINKS, replace=FALSE)
   
   sig_link_ids <- setdiff(remaining_link_ids, inh_link_ids) 
   
   # test for mutually exclusivity: 
   is12 <- intersect(act_link_ids, inh_link_ids)
   length(is12)
   
   is15 <- intersect(act_link_ids, sig_link_ids)
   length(is15)
   
   is25 <- intersect(inh_link_ids, sig_link_ids)
   length(is25)
   
   # find act type positions:
   act_pos_vec <- tfnet$LINK_ID %in% act_link_ids  
   
   # find inh type positions:
   inh_pos_vec <- tfnet$LINK_ID %in% inh_link_ids  
   
   # find sig type positions:
   sig_pos_vec <- tfnet$LINK_ID %in% sig_link_ids
   
   # assign appropriate LINK TYPE to reach link: 
   tfnet$TYPE <- replace(act_pos_vec, act_pos_vec == TRUE, TYPE_ACT) +  
      replace(inh_pos_vec, inh_pos_vec == TRUE, TYPE_INH) +
      replace(sig_pos_vec, sig_pos_vec == TRUE, TYPE_SIG)
   
   tfnet$LINK_ID <- NULL
   #return(list(tfnet, reg_db))
   return(tfnet)
}

#---------------------------------------------------------------------#
#           create TF to TARGET network
#---------------------------------------------------------------------#
create_tf2target_net <- function(TFSET, TARGETSET, 
                                 NO_TARGETS_PER_TF){
   # create a target pool and establish TF to TARGET interactions
   #NO_TARGETS <- 200 
   #NO_TARGETS_PER_TF <- 20 
   TYPE_ACT <- 1
   TYPE_INH <- 2
   
   NO_TARGETS <- length(TARGETSET) # target poolsize
   TOTAL_TF2TARGET_LINKS <- NO_TFS*NO_TARGETS_PER_TF
   
   PERCENT_ACT_LINKS_TARGET <- 0.50 # percent of activation links
   PERCENT_INH_LINKS_TARGET <- 0.50 # percent of inhibition links
   
   NO_ACT_LINKS_TARGET <- TOTAL_TF2TARGET_LINKS * PERCENT_ACT_LINKS_TARGET
   NO_INH_LINKS_TARGET <- TOTAL_TF2TARGET_LINKS * PERCENT_INH_LINKS_TARGET
   
   
   #TARGETSET <- sprintf("tg%d", 1:NO_TARGETS) 
   
   
   tf2target_net <- data.frame('SOURCE' = character(), 
                               'TARGET' = character())
   
   # for each gene, add 20 random links: 
   for (g in TFSET){ 
      targets <- sample(TARGETSET, NO_TARGETS_PER_TF)
      tf2target_subnet <- data.frame('SOURCE' = rep(g, each=NO_TARGETS_PER_TF), 
                                     'TARGET' = targets)
      tf2target_net <- rbind(tf2target_net, tf2target_subnet)
      
      # add the new members to the regulon DB:
      #reg_db[[g]] <- c(reg_db[[g]], targets)
   }
   
   # add a column LINK_ID: 
   tf2target_net <- cbind(LINK_ID = seq(1:TOTAL_TF2TARGET_LINKS),  
                          tf2target_net)
   
   #tf2target_net$LINK_ID <- NULL
   #tf2target_net$TYPE <- NULL
   
   # add a column TYPE: 
   tf2target_net <- cbind(tf2target_net, TYPE = 0)
   
   act_link_ids <- sample(tf2target_net$LINK_ID, 
                          NO_ACT_LINKS_TARGET, 
                          replace=FALSE) 
   
   inh_link_ids <- setdiff(tf2target_net$LINK_ID, 
                           act_link_ids) 
   
   
   # test for mutually exclusivity: 
   is12 <- intersect(act_link_ids, inh_link_ids)
   length(is12)
   
   
   # find act type positions:
   act_pos_vec <- tf2target_net$LINK_ID %in% act_link_ids  
   
   # find inh type positions:
   inh_pos_vec <- tf2target_net$LINK_ID %in% inh_link_ids  
   
   # assign appropriate LINK TYPE to reach link: 
   tf2target_net$TYPE <- replace(act_pos_vec, act_pos_vec == TRUE, TYPE_ACT) +  
      replace(inh_pos_vec, inh_pos_vec == TRUE, TYPE_INH)  
   
   tf2target_net$LINK_ID <- NULL   
   #return(list(tf2target_net, reg_db))
   return(tf2target_net)
}

#---------------------------------------------------------------------#
#           create regulon DB
#---------------------------------------------------------------------#
create_regdb <- function(tfnet){
   reg_db <- list()
   TYPE_ACT <- 1
   TYPE_INH <- 2
   TYPE_SIG <- 5
   
   for(i in seq(1:dim(tfnet)[1])){
      source <- as.character(tfnet[i,"SOURCE"])
      target <- as.character(tfnet[i,"TARGET"])
      type <- as.character(tfnet[i,"TYPE"])
      
      if (type == TYPE_SIG) 
         next
      if(!length(reg_db[[source]])) 
         reg_db[[source]] <- target
      else 
         reg_db[[source]] <- c(reg_db[[source]], target)
   }
   return(reg_db)
}

#---------------------------------------------------------------------#
# calculate connectivity matrix using the regulon DB
#---------------------------------------------------------------------#
create_conn_matrix <- function(reg_db){
   source_list <- names(reg_db)
   target_list <- unique(unlist(reg_db))
   
   library(gtools)
   # order the targets lexicographically: 
   source_list <- source_list[mixedorder(source_list)]
   
   # order the targets lexicographically: 
   target_list <- target_list[mixedorder(target_list)]
   
   conn_df <- data.frame(matrix(nrow = length(target_list), 
                                ncol=length(source_list)))
   conn_df[is.na(conn_df)] <- 0
   row.names(conn_df) <- target_list
   colnames(conn_df) <- source_list
   
   # place 1 when a target is found for each TF: 
   for (tf in names(reg_db)){
      for (tg in reg_db[[tf]]) {
         conn_df[tg, tf] <- 1 
      }
   }
   #return(conn_df)
   return(list(conn_df, target_list))
}


#---------------------------------------------------------------------#
# calculate correlation from ndata for specific type of interaction
#---------------------------------------------------------------------#
cal_corr <- function(ndata, tfnet, TYPE){
   cor_arr <- vector(length = 0)
   
   for (i in seq(1, dim(tfnet)[1])) { 
      type <- as.character(tfnet[i,"TYPE"])
      gene_pair <- c((as.character(tfnet$SOURCE))[i], 
                     (as.character(tfnet$TARGET))[i])
      if(type==TYPE){
         cor_arr <- c(cor_arr, cor(ndata[,gene_pair[1]], 
                                   ndata[,gene_pair[2]], 
                                   method = "spearman")) 
      }
   }
   return(cor_arr)
}

# consider all interaction types: 
cal_corr_all <- function(ndata, tfnet){
   cor_arr <- vector(length = 0)
   
   for (i in seq(1, dim(tfnet)[1])) { 
      gene_pair <- c((as.character(tfnet$SOURCE))[i], 
                     (as.character(tfnet$TARGET))[i])
         cor_arr <- c(cor_arr, cor(ndata[,gene_pair[1]], 
                                   ndata[,gene_pair[2]], 
                                   method = "spearman")) 
   }
   return(cor_arr)
}

# consider specific interaction types - ndata.states vs ndata.exp: 
cal_corr_states_vs_exp <- function(ndata.states, ndata.exp, tfnet, TYPE){
   cor_arr <- vector(length = 0)
   
   for (i in seq(1, dim(tfnet)[1])) { 
      type <- as.character(tfnet[i,"TYPE"])
      gene_pair <- c((as.character(tfnet$SOURCE))[i], 
                     (as.character(tfnet$TARGET))[i])
      if(type==TYPE){
         cor_arr <- c(cor_arr, cor(ndata.states[,gene_pair[1]], 
                                   ndata.exp[,gene_pair[2]], 
                                   method = "spearman"))
      }
   }
   return(cor_arr)
}

# consider all interaction types - ndata.states vs ndata.exp: 
cal_corr_all_states_vs_exp <- function(ndata.states, ndata.exp, tfnet){
   cor_arr <- vector(length = 0)
   
   for (i in seq(1, dim(tfnet)[1])) { 
      gene_pair <- c((as.character(tfnet$SOURCE))[i], 
                     (as.character(tfnet$TARGET))[i])
      cor_arr <- c(cor_arr, cor(ndata.states[,gene_pair[1]], 
                                ndata.exp[,gene_pair[2]], 
                                method = "spearman")) 
   }
   return(cor_arr)
}

"
   Calculate correlations between target pairs of each TF 
"
cal_cor_between_target_pairs <- function (ndata.states, ndata.exp, tfnet){
   cor_arr <- vector(length = 0) 
   tf_set <- unique(tfnet$SOURCE)
   for (i in seq(1, length(tf_set))){
      tf <- as.character(tf_set[i])
      
      tmp <- tfnet[tfnet$SOURCE==tf,]
      target_set <- as.character(tmp$TARGET)
      target_set <- setdiff(target_set, tf_set)
      
      for(j in seq(1, (length(target_set)-1))) {
         for(k in seq(j+1, length(target_set))) { 
            gene_pair <- c(target_set[j], target_set[k]) 
            cor_arr <- c(cor_arr, cor(ndata.states[,gene_pair[1]], 
                                      ndata.exp[,gene_pair[2]], 
                                      method = "spearman")) 
         }
      }   
   }   
   return(cor_arr)
}

"
   Calculate mitual information between target pairs of each TF 
"
cal_mi_between_target_pairs <- function (ddata.states, ddata.exp, tfnet){
   
   mi_arr <- vector(length = 0) 
   tf_set <- unique(tfnet$SOURCE)
   for (i in seq(1, length(tf_set))){
      tf <- as.character(tf_set[i])
      
      tmp <- tfnet[tfnet$SOURCE==tf,]
      target_set <- as.character(tmp$TARGET)
      target_set <- setdiff(target_set, tf_set)
      
      for(j in seq(1, (length(target_set)-1))) {
         for(k in seq(j+1, length(target_set))) { 
            gene_pair <- c(target_set[j], target_set[k]) 
            mi_arr <- c(mi_arr, mutinformation(ddata.states[,gene_pair[1]], 
                                               ddata.exp[,gene_pair[2]], 
                                               method = "emp")) 
         }
      }   
   }
   return(mi_arr)
}


#---------------------------------------------------------------------#
# calculate mitual information from 
# discretized data for specific type of interaction
#---------------------------------------------------------------------#
cal_mi <- function(ddata, tfnet, TYPE){
   mi_arr <- vector(length = 0)
   
   for (i in seq(1, dim(tfnet)[1])) { 
      type <- as.character(tfnet[i,"TYPE"])
      gene_pair <- c((as.character(tfnet$SOURCE))[i], 
                     (as.character(tfnet$TARGET))[i])
      if(type==TYPE){
         mi_arr <- c(mi_arr, mutinformation(ddata[,gene_pair[1]], 
                                            ddata[,gene_pair[2]], 
                                            method = "emp")) 
      }
   }
   return(mi_arr)
}

# consider all links: 
cal_mi_all <- function(ddata, tfnet){
   mi_arr <- vector(length = 0)
   for (i in seq(1, dim(tfnet)[1])) { 
      gene_pair <- c((as.character(tfnet$SOURCE))[i], 
                     (as.character(tfnet$TARGET))[i])
      mi_arr <- c(mi_arr, mutinformation(ddata[,gene_pair[1]], 
                                         ddata[,gene_pair[2]], 
                                         method = "emp")) 
   }
   return(mi_arr)
}

# calculate MI for specific type of interactions - states vs exp:
cal_mi_states_vs_exp <- function(ddata.states, ddata.exp, tfnet, TYPE){
   mi_arr <- vector(length = 0)
   
   for (i in seq(1, dim(tfnet)[1])) { 
      type <- as.character(tfnet[i,"TYPE"])
      gene_pair <- c((as.character(tfnet$SOURCE))[i], 
                     (as.character(tfnet$TARGET))[i])
      if(type==TYPE){
         mi_arr <- c(mi_arr, mutinformation(ddata.states[,gene_pair[1]], 
                                            ddata.exp[,gene_pair[2]], 
                                            method = "emp")) 
      }
   }
   return(mi_arr)
}


# consider all links - states vs exp: 
cal_mi_all_states_vs_exp <- function(ddata.states, ddata.exp, tfnet){
   mi_arr <- vector(length = 0)
   for (i in seq(1, dim(tfnet)[1])) { 
      gene_pair <- c((as.character(tfnet$SOURCE))[i], 
                     (as.character(tfnet$TARGET))[i])
      mi_arr <- c(mi_arr, mutinformation(ddata.states[,gene_pair[1]], 
                                         ddata.exp[,gene_pair[2]], 
                                         method = "emp")) 
   }
   return(mi_arr)
}

#---------------------------------------------------------------------#
#           calculate frequencies
#---------------------------------------------------------------------#
cal_frequency <- function(elements) {
   elements_sorted  <- sort(elements)
   elements_dup <- duplicated(elements_sorted)
   pos_vector <- which(elements_dup ==FALSE)
   freq_vector <- vector(length=0)
   for (i in seq(1, length(pos_vector))){
      freq_vector <- c(freq_vector, pos_vector[i]-pos_vector[i-1])
   }
   return(freq_vector)
}

#---------------------------------------------------------------------#
# plot histogram for correlations 
#---------------------------------------------------------------------#
plot_cor_hist_TF30 <- function(fig_name=fig_name, 
                               cor_random, 
                               cor_all, 
                               cor_act, 
                               cor_sig, 
                               cor_inh){
   breaks = seq(-1, 1, by=0.1)
   
   pdf(file = fig_name, width=3, height=8, 
       paper='special', onefile = TRUE)
   par(mfrow = c(5,1))
   hist(cor_random, breaks = breaks, xlim = c(-1.0, 1.0), 
        freq = FALSE, main = 'non-interactions')
   par(new=FALSE)
   hist(cor_all, breaks = breaks,
        xlim = c(-1.0, 1.0), freq = FALSE, main = 'all interactions')
   par(new=FALSE)
   hist(cor_act, breaks = breaks,
        xlim = c(-1.0, 1.0), freq = FALSE, main = 'excitatory')
   par(new=FALSE)
   hist(cor_sig, breaks = breaks,
        xlim = c(-1.0, 1.0), freq = FALSE, main = 'signaling')
   par(new=FALSE)
   hist(cor_inh, breaks = breaks,
        xlim = c(-1.0, 1.0), freq = FALSE, main = 'inhibitory')
   dev.off()
}


plot_cor_hist_TF_vs_TG <- function(fig_name=fig_name, 
                               cor_random, 
                               cor_all, 
                               cor_act, 
                               cor_sig, 
                               cor_inh,  
                               cor_tg_pairs) {
   breaks = seq(-1, 1, by=0.1)
   
   pdf(file = fig_name, width=2.5, height=5.0, 
       paper='special', onefile = TRUE)
   par(mfrow = c(5,1))
   
   par(mar = c(5, 5, 0, 0), oma = c(1, 1, 1, 1))
   #par(mar = c(5, 5, 0, 0), oma = c(0, 0, 0, 0))
   par(mgp = c(2.5, 1.0, 0))
   par(tcl = -0.25)
   par(cex = 0.4, cex.lab=2.0, cex.main=2.0) 

   xlim <- c(-1.0, 1.0)
   ylim <- c(0.0, 3.2)
   at.yaxis = c(0.0, 1.5, 3.0)
   xtic.font.size = 1.8
   ytic.font.size = 1.8
   fil.color = 'gray'
   hist(cor_random, breaks = breaks, col = fil.color,
        xlab = "",  ylab = "",
        xaxt='n', yaxt='n',
        xlim = xlim, ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at=at.yaxis) 
   
   par(new=FALSE)
   hist(cor_all, breaks = breaks, col = fil.color,
        xlab = "", ylab = "", 
        xaxt='n', yaxt='n',
        xlim = xlim, ylim = ylim, 
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at=at.yaxis) 
   
   par(new=FALSE)
   hist(cor_act, breaks = breaks, col = fil.color,
        xlab = "", ylab = "",  
        xaxt='n', yaxt='n',
        xlim = xlim, ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   
   par(new=FALSE)
   hist(cor_sig, breaks = breaks, col = fil.color,
        xlab = "",  ylab = "", 
        xaxt='n',  yaxt='n',
        xlim = xlim, ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   
   par(new=FALSE)
   hist(cor_inh, breaks = breaks, col = fil.color, 
        xlab = "", ylab = "", 
        xaxt='n',  yaxt='n',
        xlim = xlim, ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   axis(1, cex.axis=xtic.font.size, at = c(-1.0, -0.5, 0, 0.5, 1.0), 
        labels = TRUE)
   dev.off()
}

plot_cor_hist_TF_vs_TG.kdtf9 <- function(fig_name=fig_name, 
                                   cor_random, 
                                   cor_all, 
                                   cor_act, 
                                   cor_sig, 
                                   cor_inh,  
                                   cor_tg_pairs) {
   breaks = seq(-1, 1, by=0.1)
   
   pdf(file = fig_name, width=2.3, height=5.0, 
       paper='special', onefile = TRUE)
   par(mfrow = c(5,1))
   
   par(mar = c(5, 5, 0, 0), oma = c(1, 1, 1, 1))
   #par(mar = c(5, 5, 0, 0), oma = c(0, 0, 0, 0))
   #par(mgp = c(2.5, 0.6, 0))
   par(mgp = c(2.5, 1.5, 0))
   
   par(tcl = -0.25)
   par(cex = 0.4, cex.lab=2.0, cex.main=2.0) 
   
   xlim <- c(-1.0, 1.0)
   ylim <- c(0.0, 3.2)
   xtic.font.size = 1.8
   ytic.font.size = 1.8
   
   at.yaxis = c(0.0, 1.5, 3.0)
   at.xaxis = c(-1.0, -0.5, 0, 0.5, 1.0)
   
   fil.color = 'gray'
   hist(cor_random, breaks = breaks, col = fil.color,
        xlab = "",  ylab = "",
        xaxt='n', yaxt='n',
        xlim = xlim, ylim = ylim,
        freq = FALSE, main = "")
   #axis(2, cex.axis=ytic.font.size) 
   axis(2, cex.axis=ytic.font.size, at=at.yaxis)
   
   par(new=FALSE)
   hist(cor_all, breaks = breaks, col = fil.color,
        xlab = "", ylab = "", 
        xaxt='n', yaxt='n',
        xlim = xlim, ylim = ylim, 
        freq = FALSE, main = "")
   #axis(2, cex.axis=ytic.font.size) 
   axis(2, cex.axis=ytic.font.size, at=at.yaxis) 
   
   par(new=FALSE)
   hist(cor_act, breaks = breaks, col = fil.color,
        xlab = "", ylab = "",  
        xaxt='n', yaxt='n',
        xlim = xlim, ylim = ylim,
        freq = FALSE, main = "")
   #axis(2, cex.axis=ytic.font.size) 
   axis(2, cex.axis=ytic.font.size, at=at.yaxis) 
   
   par(new=FALSE)
   hist(cor_sig, breaks = breaks, col = fil.color,
        xlab = "",  ylab = "", 
        xaxt='n',  yaxt='n',
        xlim = xlim, ylim = ylim,
        freq = FALSE, main = "")
   #axis(2, cex.axis=ytic.font.size) 
   axis(2, cex.axis=ytic.font.size, at=at.yaxis) 
   
   par(new=FALSE)
   hist(cor_inh, breaks = breaks, col = fil.color, 
        xlab = "", ylab = "", 
        xaxt='n',  yaxt='n',
        xlim = xlim, ylim = ylim,
        freq = FALSE, main = "")
   #axis(2, cex.axis=ytic.font.size) 
   axis(2, cex.axis=ytic.font.size, at=at.yaxis) 
   axis(1, cex.axis=xtic.font.size, at = at.xaxis, 
        labels = TRUE)
   dev.off()
}




plot_cor_hist_targetpairs <- function(fig_name=fig_name, 
                                   cor_tg_pairs) {
   breaks = seq(-1, 1, by=0.1)
   
   pdf(file = fig_name, width=2.5, height=1.2, 
       paper='special', onefile = TRUE)
   par(mfrow = c(1,1))
   
   par(mar = c(5, 5, 0, 0), oma = c(1, 1, 1, 1))
   #par(mar = c(5, 5, 0, 0), oma = c(0, 0, 0, 0))
   par(mgp = c(2.5, 1.0, 0))
   par(tcl = -0.25)
   par(cex = 0.4, cex.lab=2.0, cex.main=2.0) 
   
   xlim <- c(-1.0, 1.0)
   ylim <- c(0.0, 3.0)
   at.yaxis = c(0.0, 1.5, 3.0)
   xtic.font.size = 1.8
   ytic.font.size = 1.8
   fil.color = 'gray'

   hist(cor_tg_pairs, breaks = breaks, col = fil.color, 
        xlab = "", ylab = "", 
        xaxt='n',  yaxt='n',
        xlim = xlim, ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   axis(1, cex.axis=xtic.font.size, at = c(-1.0, -0.5, 0, 0.5, 1.0), 
        labels = TRUE)
   dev.off()
}









# plot histogram for mitual information
plot_mi_hist_TF30 <- function(fig_name=fig_name, mi_random, 
                         mi_all, mi_act, mi_sig, mi_inh) {
   breaks = seq(0, 1.44, by=0.03)
   xlim = c(0, 0.5)
   
   pdf(file = fig_name, width=3, height=8, 
       paper='special', onefile = TRUE)
   
   
   par(mfrow = c(5,1))
   hist(mi_random, breaks = breaks, xlim = xlim,  
        freq = FALSE, main = "non-interactions")
   par(new=FALSE)
   
   hist(mi_all, breaks = breaks, xlim = xlim, 
        freq = FALSE, main = "all interactions")
   par(new=FALSE)
   
   hist(mi_act, breaks = breaks, xlim = xlim, 
        freq = FALSE, main = "excitatory")
   par(new=FALSE)
   
   hist(mi_sig, breaks = breaks, xlim = xlim, 
        freq = FALSE, main = "signaling")
   par(new=FALSE)
   
   hist(mi_inh, breaks = breaks, xlim = xlim, 
        freq = FALSE, main = "inhibitory")
   dev.off()
}


plot_mi_hist_TF30.nice <- function(fig_name=fig_name, mi_random, 
                              mi_all, mi_act, mi_sig, mi_inh) {
  
   
   pdf(file = fig_name, width=2.5, height=5.0, 
       paper='special', onefile = TRUE)
   
   par(mfrow = c(5,1))
   
   par(mar = c(5, 5, 0, 0), oma = c(1, 1, 1, 1))
   #par(mgp = c(2.5, 0.6, 0))
   par(mgp = c(2.5, 1.5, 0))
   #par(mgp = c(0, 0.6, 0))
   
   par(tcl = -0.25)
   par(cex = 0.4, cex.lab=2.0, cex.main=2.0) 
   
   xlim = c(0, 0.6)  
   ylim = c(0, 16.0)
   xtic.font.size = 1.8
   ytic.font.size = 1.8
   fil.color = 'gray'
   
   at.yaxis = c(0, 8, 16)  
   at.xaxis = c(0,  0.2,  0.4, 0.6)
   
   breaks = seq(0, 1.44, by=0.03)

   # main = "non-interactions"
   hist(mi_random, breaks = breaks, 
        xlab = "",  ylab = "",
        xaxt='n',  yaxt='n',
        col = fil.color,
        xlim = xlim,  ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   
   par(new=FALSE)
   # main = "all interactions"
   hist(mi_all, breaks = breaks, 
        xaxt='n',  yaxt='n',
        xlab = "",  ylab = "",
        col = fil.color,
        xlim = xlim,  ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   
   par(new=FALSE)
   # main = "excitatory"
   hist(mi_act, breaks = breaks, 
        xaxt='n',  yaxt='n',
        xlab = "",  ylab = "", 
        col = fil.color,
        xlim = xlim,  ylim = ylim,
        freq = FALSE, main = "" )
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   
   par(new=FALSE)
   # main = "signaling"
   hist(mi_sig, breaks = breaks, 
        xaxt='n',  yaxt='n',
        xlab = "",  ylab = "",
        col = fil.color,
        xlim = xlim,  ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   
   par(new=FALSE)
   # main = "inhibitory"
   hist(mi_inh, breaks = breaks, 
        xaxt='n',  yaxt='n',     
        xlab = "",  ylab = "",
        col = fil.color,
        xlim = xlim,  ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   axis(1, cex.axis=xtic.font.size, at = at.xaxis, 
        labels = TRUE)   
   
   dev.off()
}


plot_mi_hist_TF30.nice.kdtf9 <- function(fig_name=fig_name, mi_random, 
                                   mi_all, mi_act, mi_sig, mi_inh) {
   
   pdf(file = fig_name, width=2.3, height=5.0, 
       paper='special', onefile = TRUE)
   
   par(mfrow = c(5,1))
   
   par(mar = c(5, 5, 0, 0), oma = c(1, 1, 1, 1))
   #par(mgp = c(2.5, 0.6, 0))
   #par(mgp = c(2.5, 1.4, 0))
   par(mgp = c(2.5, 1.5, 0))
   
   #par(mgp = c(0, 0.6, 0))
   
   par(tcl = -0.25)
   par(cex = 0.4, cex.lab=2.0, cex.main=2.0) 
   
   #xlim = c(0, 0.5)
   xlim = c(0, 0.65)
   ylim = c(0, 9.0)
   
   xtic.font.size = 1.8
   ytic.font.size = 1.8
   fil.color = 'gray'
   
   at.yaxis = c(0,  4,  8) 
   at.xaxis = c(0,  0.2,  0.4,  0.6)
   
   breaks = seq(0, 1.44, by=0.03)
   
   # main = "non-interactions"
   hist(mi_random, breaks = breaks, 
        xlab = "",  ylab = "",
        xaxt='n',  yaxt='n',
        col = fil.color,
        xlim = xlim,  ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   
   par(new=FALSE)
   # main = "all interactions"
   hist(mi_all, breaks = breaks, 
        xaxt='n',  yaxt='n',
        xlab = "",  ylab = "",
        col = fil.color,
        xlim = xlim,  ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   
   par(new=FALSE)
   # main = "excitatory"
   hist(mi_act, breaks = breaks, 
        xaxt='n',  yaxt='n',
        xlab = "",  ylab = "",   
        col = fil.color,
        xlim = xlim,  ylim = ylim,
        freq = FALSE, main = "" )
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   
   par(new=FALSE)
   # main = "signaling"
   hist(mi_sig, breaks = breaks, 
        xaxt='n',  yaxt='n',
        xlab = "",  ylab = "",
        col = fil.color,
        xlim = xlim,  ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   
   par(new=FALSE)
   # main = "inhibitory"
   hist(mi_inh, breaks = breaks, 
        xaxt='n',  yaxt='n',     
        xlab = "",  ylab = "",
        col = fil.color,
        xlim = xlim,  ylim = ylim,
        freq = FALSE, main = "")
   axis(2, cex.axis=ytic.font.size, at = at.yaxis) 
   axis(1, cex.axis=xtic.font.size, at = at.xaxis, 
        labels = TRUE)   
   
   dev.off()
}





#---------------------------------------------------------------------#
#  Extract activaty from RACIPE and FastNCA
#---------------------------------------------------------------------#
xTract_activities <- function(act_nca,  act_racipe){
   act_pairs <- list()
   tf_list <- colnames(act_nca)

   for (tf in tf_list){      
      # extract RACIPE activities: 
      v1 <- as.numeric(act_racipe[, tf])
      
      # extract FastNCA activities:
      v2 <- as.numeric(act_nca[, tf])
      
      # insert two activities as a pair in the list:
      act_pair <- list('racipe'=v1, 
                       'nca'=v2)
      act_pairs[[tf]] <-  act_pair
   }
   return(act_pairs)
}

xTract_activities_old <- function(se_df,  conn_df){
   act_pairs <- list()
   for (i in seq(1, dim(se_df)[1])){
      tf_name <- colnames(conn_df)[i] 
      
      # extract RACIPE activities: 
      v1 <- act_df[, tf_name]
      
      # extract FastNCA activities:
      v2 <- as.numeric(se_df[i,])
      
      # insert two activities as a pair in the list:
      act_pair <- list('racipe'=v1, 
                       'nca'=v2)
      act_pairs[[tf_name]] <-  act_pair
   }
   return(act_pairs)
}

#---------------------------------------------------------------------#
# Calculate correlations from activation pairs from RACIPE and FastNCA
#---------------------------------------------------------------------#
cal_correlations <- function(act_pairs){
   cor_vector <- vector(length=0)
   tf_list <- names(act_pairs)
   
   for (tf in tf_list){    
      v1 <- act_pairs[[tf]][['racipe']]
      v2 <- act_pairs[[tf]][['nca']]

      cor_vector <- c(cor_vector, 
                      cor(v1, v2, method = 'spearman'))
   }
   return(cor_vector)
}

#---------------------------------------------------------------------#
# Plor activation pairs from RACIPE and FastNCA
#---------------------------------------------------------------------#
plot_activities_racipe_vs_nca <-function(act_pairs, fname_fig){
   NO_WINDOWS_PER_PAGE <-5
   
   pdf(file = fname_fig, width=10, height=20, 
       paper='special', onefile = TRUE)
   tf_list <- names(act_pairs)  
   
   par(mfrow=c(NO_WINDOWS_PER_PAGE,1))
   
   for (tf in tf_list){
      v1 <- act_pairs[[tf]][['racipe']]
      v2 <- act_pairs[[tf]][['nca']]

      v1 <- (v1-min(v1))/(max(v1) - min(v1))
      v2 <- (v2-min(v2))/(max(v2) - min(v2))

      plot(v1 , col='black', xaxt = "n",  
           xlab = "", ylab = "", type = 'b', cex.lab=2.5)
      par(new=TRUE)
      plot(v2 , col='red', xlab = tf, 
           ylab = "", type = 'b', cex.lab=2.5)
      #break
   }
   dev.off()
}


plot_activities_racipe_vs_nca_old <-function(act_pairs, fname_fig){
   pdf(file = fname_fig, width=10, height=20, 
       paper='special', onefile = TRUE)
   par(mfrow=c(dim(se_df)[1]/5,1))
   tf_list <- names(act_pairs)
   for (tf in tf_list){
      v1 <- act_pairs[[tf]][['racipe']]
      v2 <- act_pairs[[tf]][['nca']]
      
      v1 <- (v1-min(v1))/(max(v1) - min(v1))
      v2 <- (v2-min(v2))/(max(v2) - min(v2))
      
      plot(v1 , col='black', xaxt = "n",  
           xlab = "", ylab = "", type = 'b', cex.lab=2.5)
      par(new=TRUE)
      plot(v2 , col='red', xlab = tf, 
           ylab = "", type = 'b', cex.lab=2.5)
   }
   dev.off()
}

#======================= for NCA correlation distributions =============================

#---------------------------------------------------------------------#
# Create a perturbed regulon DB for a supplied one
#---------------------------------------------------------------------#
create_perturbed_regdb <- function(regdb.no, 
                                   NO_TARGETES_2B_PERTURBED){
   tflist <- names(regdb.no)
   targets_all <- unique(unlist(regdb.no))
   
   target_pool <- setdiff(targets_all, tflist)
   
   regdb <- list()
   
   for (tf in tflist){
      targets <- regdb.no[[tf]]
      target_pool_cur <- setdiff(target_pool,  targets)
      
      # TFs in the targets for the current TF:
      tfs_in_targets <- setdiff(targets, target_pool)
      
      targets_remaining <- setdiff(targets, tfs_in_targets)
      
      # randomly select targets to be removed
      targets_2b_removed <- sample(targets_remaining, 
                                   NO_TARGETES_2B_PERTURBED) 
      
      # remove selected targets: 
      targets_remaining <- setdiff(targets_remaining, 
                                   targets_2b_removed)
      # ranomly select targets to be added from the target pool for the current TF: 
      targets_2b_added <- sample(target_pool_cur, 
                                 NO_TARGETES_2B_PERTURBED)
      # create targets for the current TF: 
      targets_new <- c(tfs_in_targets, targets_remaining, targets_2b_added)
      
      # insert into the regulon DB:
      regdb[[tf]] <- targets_new
   }   
   return(regdb)
}


#---------------------------------------------------------------------#
#  Histogram for correaltions across perturbed regulon DBs of specific size
#---------------------------------------------------------------------#
plot_cor.across_regDBs <- function(cor_df, fname_fig) {
   #m <- t(as.matrix(cor_df)) 
   m <- as.matrix(cor_df)
   cor_all <- as.vector(m)
   
   cor_rowMeans <- rowMeans(abs(m))
   cor_colMeans <- colMeans(abs(m))
   
   NO_WINDOWS_PER_PAGE <- 3
   
   pdf(file = fname_fig, width=4, height=8, 
       paper='special', onefile = TRUE)
   par(mfrow=c(NO_WINDOWS_PER_PAGE,1))
   
   breaks = seq(0, 1.0, by=0.01)
   hist(abs(cor_all), breaks = breaks)
   
   breaks = seq(0, 1.0, by=0.01)
   hist(cor_rowMeans, breaks = breaks, main = 'row means: mean across TFs')
   
   breaks = seq(0, 1.0, by=0.01)
   hist(cor_colMeans, breaks = breaks, main = 'col means: mean across regulon DBs')
   dev.off()
}

