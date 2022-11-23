"
   Creates the following files: 
      (1) net30tf.states.txt 
         It creates this activation file by combining 
         the activation files from racipe simulations of 
         wild and mutant type network.
      (2) net30tf.states.exp.txt
         It creates this expression file by combining 
         the activation files from racipe simulations of 
         wild and mutant type network.
"

#   Function definition
#======================#
"
   Creates the following files by combining wild and mutant type data: 
      (1) net30tf.states.txt 
         It creates this activation file by combining 
         the activation files from racipe simulations of 
         wild and mutant type network.
      (2) net30tf.states.exp.txt
         It creates this expression file by combining 
         the activation files from racipe simulations of 
         wild and mutant type network.
"

comb.racipe.states.mt_and_wt <- function(){ 
  #-------------------------------------------------------#
  #   CONSTANTS
  #-------------------------------------------------------#
  # number of models to be selected from wild type simulation:
  MODELS_WT <- 40 
  # number of models to be selected from mutant type simulation:
  MODELS_MT <- 43
  
  #-----------------------------------------------#
  #   Define directory variables
  #-----------------------------------------------#
  datadir_wt <- '../sim.racipe.wt/'
  datadir_mt <- '../sim.racipe.kd9/'
  
  #-----------------------------------------------#
  #   Define output directories
  #-----------------------------------------------#
  combdir <- "./"
  #dir.create(combdir)
  
  #------------------------------------------------#
  #   Define input filenames
  #------------------------------------------------#
  # regulon DB file name: 
  fname_regulon <- paste0(datadir_wt, 'regdb.net30tf.perturbed10.rda')  
  
  # Activities files
  fname_act_wt <- paste0(datadir_wt, 'net30tf.states.txt') 
  fname_act_mt <- paste0(datadir_mt, 'net30tf.states.txt') 
  
  # Expression files
  fname_exp_wt <- paste0(datadir_wt, 'net30tf.states.exp.txt') 
  fname_exp_mt <- paste0(datadir_mt, 'net30tf.states.exp.txt') 
  
  #------------------------------------------------------#
  #   Define output filenames - for FastNCA input format
  #------------------------------------------------------#
  # combined activity file:
  fname_act_comb_txt <- paste0(combdir, 'net30tf.states.txt')
  
  # combined expression file:
  fname_exp_comb_txt <- paste0(combdir, 'net30tf.states.exp.txt')
  
  #--------------------------------------------------------#
  #   Load RACIPE data - activities and expressions
  #--------------------------------------------------------#
  # activity files:
  act_wt <- read.table(fname_act_wt, header = TRUE)
  act_mt <- read.table(fname_act_mt, header = TRUE)
  
  # expression files:
  exp_wt <- read.table(fname_exp_wt, header = TRUE)
  exp_mt <- read.table(fname_exp_mt, header = TRUE)
  
  # create row names:
  obs_ids <- c(paste0('M', act_wt$MODEL_NO[1:MODELS_WT]), 
               paste0('T', act_mt$MODEL_NO[1:MODELS_MT]))
  
  # combine activities from wild and mutant types:
  act_comb <- rbind(act_wt[1:MODELS_WT, ], 
                    act_mt[1:MODELS_MT,])
  act_comb$MODEL_NO <- obs_ids
  
  # combine expressions from wild and mutant types:
  exp_comb <- rbind(exp_wt[1:MODELS_WT, ], exp_mt[1:MODELS_MT,])
  exp_comb$MODEL_NO <- obs_ids
  
  #---------------------------------------------------------#
  #   Save combined activities in RACIPE Format: as a text and as a csv file
  #---------------------------------------------------------#
  write.table(act_comb, file = fname_act_comb_txt, 
              quote = FALSE, row.names = FALSE)
  
  write.table(exp_comb, file = fname_exp_comb_txt, 
              quote = FALSE, row.names = FALSE)
  return(0)
}

comb.racipe.states.mt_and_wt()
