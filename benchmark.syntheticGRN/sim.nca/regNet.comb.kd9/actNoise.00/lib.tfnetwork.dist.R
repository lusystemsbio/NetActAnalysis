#---------------------------------------------------------------------#
# Calculate connectivity matrix using the regulon DB
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
      for (tg in reg_db[[tf]]){
         conn_df[tg, tf] <- 1 
      }
   }
   
   conn.mat <- list()
   conn.mat[["CONN.MAT"]] <- conn_df
   conn.mat[["TARGET.LIST"]] <- target_list

   return(conn.mat) 
}


cal.nca.act <- function(regdb.list, exp_df){
   cat("\nCalculating NCA activities for all perturbed regulon DBs ...\n") 
   nca_act.list <- list()
   #NO_REGDB_2BUSED <- 2 # this is for debug purpose
   #regdb.list <- regdb.list[1:NO_REGDB_2BUSED] # comment out for using all DBs
   
   for(regdb_name in names(regdb.list)){
      create.inputs.forFastNCA(regdb=regdb.list[[regdb_name]], exp_df=exp_df) 
      system("octave calFastNCA.m") # Calculate NCA activities 
      
      # format and import generated NCA activities:
      nca_act.list[[regdb_name]] <- format.NCAoutput()
   }
   return(nca_act.list)
}


# Create.inputs.forFastNCA(regdb)
#-----------------------------------------------#
"
   Creates input files for FastNCA simulation. 
   uses the following files: 
      (1) regulon DB - supplied as a parameter
      (2) racipe activation file - loads from racipe simulation folder 
      (3) racipe expression file - loads from racipe simulation folder 
  outputs the following files: 
      (1) Expression file in FastNCA input file format - 
               saved in FastNCA simulation folder 
      (2) Connectivity matrix in FastNCA input file format - constructed from the regulon DB - 
               saved in FastNCA simulation folder 
"

create.inputs.forFastNCA <- function(regdb, exp_df) { 
   library(gtools) 
   # Create connectivity matrix from regulon DB
   ret_list <- create_conn_matrix(regdb)
   conn_df <- ret_list[["CONN.MAT"]]
   target_list <- ret_list[["TARGET.LIST"]]
   write.csv(conn_df, file='net30tf.matrix.csv', quote = FALSE)
   
   # Change RACIPE expressions to FastNCA input format
   exp_df <- exp_df[target_list, ] # Retain the rows for targets only 
   write.csv(exp_df, file='net30tf.states.exp.csv', quote = FALSE)
   
   # # Construct return list
   # nca.inputs <- list() 
   # nca.inputs[["CONN.MAT"]] <- conn_df 
   # nca.inputs[["EXPR.MAT"]] <- exp_df
   # return(nca.inputs) 
   return()
}


# format.NCAoutput()
#------------------------#
"
   Creates formatted and annotated activities from FastNCA.m
   Input files: 
         (1) Se.csv
         (2) net30tf.states.exp.comb.csv
         (3) net30tf.matrix.csv

   Output file(s): 
         (1) nca.act.txt, nca.act.csv
"
format.NCAoutput <- function() {
   # load network tolopogy: to get the order of TF list in Se:
   conn_df <- read.csv(file='net30tf.matrix.csv', header=TRUE, row.names = 1)
   # load RACIPE expressions:
   exp_df <- read.csv(file = 'net30tf.states.exp.csv', header=TRUE, row.names = 1)  
   # load calculated NCA TFA matrix - Se file: FastNCA activities file:
   se_df <- read.csv(file = 'Se.csv', header=FALSE)
   
   # Format NCA activities output
   #-----------------------------
   model_list <- colnames(exp_df) # retrieve the ordered models as used by NCA
   TF_list <- colnames(conn_df) # retrieve the ordered TFs as used by NCA
   #nca_act <- t(se_df) # format NCA output activities
   # row.names(nca_act) <- model_list   
   # colnames(nca_act) <- TF_list
   
   # add models as the row names, TF list as the column names:
   nca_act <- se_df
   row.names(nca_act) <- TF_list 
   colnames(nca_act) <- model_list 
   
   write.table(nca_act, file = 'nca.act.txt', quote = FALSE, row.names = FALSE)
   write.csv(nca_act, file='nca.act.csv', quote = FALSE)
   return(nca_act)
} 