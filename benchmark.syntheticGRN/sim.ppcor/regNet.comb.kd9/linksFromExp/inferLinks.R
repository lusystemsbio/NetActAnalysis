library(ppcor)

outdir <- './data/'

# Load activities data
mydata <- read.table('./data/net30tf.states.tsv', sep = '\t', header = T, row.names = 1) 
inputExpr <- as.matrix(mydata)

geneNames <- rownames(inputExpr)
rownames(inputExpr) <- c(geneNames)

# Run pcor using spearman's correlation as mentioned in the PNI paper 
# Link to paper: https://www.pnas.org/content/114/23/5822

tmp.df =  pcor(x= t(inputExpr), method = "spearman")
names(tmp.df)
tmp.estimate <- tmp.df$estimate

pcorResults =  pcor(x= t(as.matrix(inputExpr)), method = "spearman")

# Write output to a file
# https://stackoverflow.com/questions/38664241/ranking-and-counting-matrix-elements-in-r
DF = data.frame(Gene1 = geneNames[c(row(pcorResults$estimate))], Gene2 = geneNames[c(col(pcorResults$estimate))]
                , corVal = c(pcorResults$estimate), pValue =  c(pcorResults$p.value))
outDF <- DF[order(DF$corVal, decreasing=TRUE), ]
outDF <- outDF[!(outDF$Gene1==outDF$Gene2), c(1, 2, 3)] # remove self ineractions

write.table(outDF, paste0(outdir, 'inferredRel.tsv'), sep = "\t", quote = FALSE, row.names = FALSE)
