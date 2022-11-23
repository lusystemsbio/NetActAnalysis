rm(list = ls()) 

miTSH <- 0.0
NBINS <-  3

source('./functions.R')
library(NetAct)

input <- outdir <- './data/' 

# Extract regulon db with no perturbation
regdb <- readRDS('../../../../../regDBs/regDB.nopert.rds')
names(regdb)

# Load activities data
mydata <- read.table('./data/net30tf.states.tsv', sep = '\t', header = T, row.names=1) 

# Infer TF-target interactions
#-----------------------------
inferredRel.df <- inferLinks(xdata=mydata, regdb=regdb, miTh = miTSH, 
                               nbins = NBINS, method = "spearman", DPI = FALSE)  
  
write.table(format(inferredRel.df, digits = 4), file = paste(outdir, 'inferredLinks.tsv', sep = ''), 
            row.names=TRUE, sep="\t", quote = F, col.names = T)  
