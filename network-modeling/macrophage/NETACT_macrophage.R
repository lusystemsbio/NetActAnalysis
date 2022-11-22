########### apply sign switching to macrophage data ##############

rm(list = ls())
# Using saved datasets - macrophageData.RDS, counts.rds, info.txt, 2019-03-07_circuit_file11c532b7432da.RDS
library("AnnotationDbi")
library("org.Mm.eg.db")
library(NetAct)
library(reshape2)
library(gplots)
library(MASS)
library(RColorBrewer)
require(circlize)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(edgeR)
library(limma)
library(sRACIPE)
library(sva)

setMethod(f="sracipeNormalize",
          signature="RacipeSE",
          definition=function(.object)
          {
            metadataTmp <- metadata(.object)
            assayDataTmp <- assays(.object)
            geneExpression <- assayDataTmp[[1]]
            geneExpression <- log2(geneExpression)
            means <- rowMeans(geneExpression)
            sds <-  apply(geneExpression, 1, sd)
            geneExpression <- sweep(geneExpression, 1, means, FUN = "-")
            geneExpression <- sweep(geneExpression, 1, sds, FUN = "/")
            assayDataTmp2 <- SimpleList(deterministic = geneExpression)
            #  geneExpression <- data.frame(geneExpression)
            #  assayDataTmp2 <- list(deterministic = geneExpression)
            tsSims <- 0
            if(!is.null(metadata(.object)$tsSimulations)){
              tsSims <- length(metadata(.object)$tsSimulations)
              tsSimulations <- assayDataTmp[2:(tsSims + 1)]
              tsSimulations <- lapply(tsSimulations,function(x) (1+x))
              tsSimulations <- lapply(tsSimulations,log2)
              #   stochasticSimulations <-
              #      lapply(stochasticSimulations, function(x) x[,is.finite(colMeans(x))])
              tsSimulations <- lapply(tsSimulations,
                                      function(x) sweep(x, 1, means, FUN = "-"))
              tsSimulations <- lapply(tsSimulations,
                                      function(x) sweep(x, 1, sds, FUN = "/"))
              assayDataTmp2 <- c(assayDataTmp2, tsSimulations)
            }
            stochSims <- 0
            if(!is.null(metadata(.object)$stochasticSimulations)){
              stochSims <- length(metadata(.object)$stochasticSimulations)
              stochasticSimulations <- assayDataTmp[2:(stochSims + 1)]
              stochasticSimulations <- lapply(stochasticSimulations,function(x) (1+x))
              stochasticSimulations <- lapply(stochasticSimulations,log2)
              #   stochasticSimulations <-
              #      lapply(stochasticSimulations, function(x) x[,is.finite(colMeans(x))])
              stochasticSimulations <- lapply(stochasticSimulations,
                                              function(x) sweep(x, 1, means, FUN = "-"))
              stochasticSimulations <- lapply(stochasticSimulations,
                                              function(x) sweep(x, 1, sds, FUN = "/"))
              assayDataTmp2 <- c(assayDataTmp2, stochasticSimulations)
            }
            if(!is.null(metadata(.object)$knockOutSimulations)){
              knockOutSimulations <- .object$knockOutSimulations
              koSims <- length(metadata(.object)$knockOutSimulations)
              knockOutSimulations <- assayDataTmp[(2+stochSims):(koSims +1 +stochSims)]
              for(ko in seq_len(length(knockOutSimulations))){
                simData <- knockOutSimulations[[ko]]
                tmpGene <- names(knockOutSimulations[ko])
                tmpMeans <- means
                tmpMeans[which(names(simData) == tmpGene)] <- 0
                tmpSds <- sds
                tmpSds[which(names(simData) == tmpGene)] <- 1
                simData <- log2(simData)
                simData[,which(names(simData) == tmpGene)] <- 0
                simData <- sweep(simData, 1, tmpMeans, FUN = "-")
                simData <- sweep(simData, 1, tmpSds, FUN = "/")
                knockOutSimulations[[ko]] <- simData
              }
              assayDataTmp2 <- c(assayDataTmp2, knockOutSimulations)
            }
            metadataTmp$normalized <- TRUE
            assays(.object) <- assayDataTmp2
            metadata(.object) <- metadataTmp
            return(.object)
          }
)

############### processing ######################
data_macro = readRDS("counts.rds")
info = read.table("info.txt", header = T)
#data_macro <- data_macro[rowSums(data_macro > 10),]
colnames(data_macro) = info$sample
grouping <- as.character(info$sample)
grouping <- substr(grouping,1,nchar(grouping)-3)
info$condition <- grouping
grouping <- as.factor(info$condition)

#grouping = info$condition
batch = as.factor(info$batch)
info$ID = NULL


t = sapply(rownames(data_macro), function(x) strsplit(x, "[.]")[[1]], 
           USE.NAMES = FALSE)
data_macro$gene_id = t[1,]
data_macro <- data_macro[!duplicated(data_macro$gene_id),]

counts <- data_macro[,1:21]
counts <- apply(counts, 2, as.integer)
rownames(counts) <- data_macro$gene_id
# convert symbols


geneid = data_macro$gene_id
genes <- AnnotationDbi::select(org.Mm.eg.db, keys=geneid, columns="SYMBOL", keytype="ENSEMBL")
genes = genes[!duplicated(genes$ENSEMBL),]

# filtering
count_data = as.data.frame(counts)
count_data$SYMBOL = genes$SYMBOL
count_data <- subset(count_data, !is.na(SYMBOL) 
                     & !duplicated(SYMBOL) 
                     & rowSums(counts > 10) >=21)

# generate count + phenotype data
rownames(count_data) = count_data$SYMBOL
count_data$SYMBOL = NULL
count_data <- as.matrix(count_data)

phenoData = new("AnnotatedDataFrame", data = data.frame(celltype = grouping,
                                                        batch = batch))
rownames(phenoData) = colnames(count_data)

compare_list = c("UT-IL4_2h", "UT-IFNg_2h", "UT-IFNg_IL4_2h","UT-IL4_4h", "UT-IFNg_4h", "UT-IFNg_IL4_4h")

DErslt = RNAseqDegs_limma(counts = count_data, phenodata = phenoData, complist = compare_list, qval = 0.05)

e = DErslt$Overall$e
rownames(phenoData) = colnames(e)

# design = model.matrix(~1, data = pData(phenoData))

design <- model.matrix(~grouping)

combat_edata = ComBat(dat = as.matrix(e), batch = batch, mod = design, par.prior=T)

neweset = ExpressionSet(assayData = as.matrix(combat_edata), 
                        phenoData = phenoData)


macrophageData = list(neweset = neweset, DErslt = DErslt)
saveRDS(macrophageData, file = "macrophageData.RDS")

##### plot PCA results ########

##########figure1
pca = prcomp(t(combat_edata), scale = TRUE)
plotData <- as.data.frame(pca$x[,1:2])
plotData$name <- rownames(pData(phenoData))
plotData$type <- c( "UT", "IL4","IFNg", "IFNg_IL4", "IL4","IFNg", "IFNg_IL4", "UT", "IL4",      "IFNg" ,    "IFNg_IL4", "IL4"    ,  "IFNg" ,    "IFNg_IL4", "UT", "IL4",  "IFNg",    "IFNg_IL4", "IL4" ,     "IFNg" ,    "IFNg_IL4")
plotData$color <- c("black", "coral", "deepskyblue", "green1","coral4", "deepskyblue3", "green4", "black", "coral", "deepskyblue", "green1","coral4", "deepskyblue3", "green4","black", "coral", "deepskyblue", "green1","coral4", "deepskyblue3", "green4" )
plotData$ftype <- factor(plotData$type)
plotData$type <- factor (substr(plotData$name,1,nchar(plotData$name)-3) ) 
plotData$hour <- (substr(plotData$name, nchar(plotData$name)-4,nchar(plotData$name)-4))
plotData$hour[plotData$hour == "U"] <- 1
plotData$batch <- factor(substr(plotData$name, nchar(plotData$name)-1,nchar(plotData$name)))

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                              
p <- ggplot(plotData, aes(x=-PC2, y=PC1, shape = plotData$batch, color = plotData$type)) +
  geom_point(size=5) + 
  labs(x="PC2 (22.6%)", y="PC1 (37.9%)", col="Treatment", shape="Batch") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p


plot.me <- plotData
plot.me$type <- as.character(plot.me$type)
plot.me$type[plot.me$type == "IL4_2h" |plot.me$type == "IL4_4h"] <- "IL4"
plot.me$type[plot.me$type == "IFNg_2h" |plot.me$type == "IFNg_4h"] <- "IFNg"
plot.me$type[plot.me$type == "IFNg_IL4_2h" |plot.me$type == "IFNg_IL4_4h"] <- "IFNg_IL4"
plot.me$type <- as.factor(plot.me$type)


ggplot(plot.me, aes(x=-PC2, y=PC1, shape = hour, color = type)) +
  geom_point(size=5) + 
  labs(x="PC2 (22.6%)", y="PC1 (37.9%)", col="Treatment", shape="Hour") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_colour_manual(values=cbbPalette)

ggplot(plot.me, aes(x=-PC2, y=PC1, shape = hour, color = type)) +
  geom_point(size=5) + 
  labs(x="PC2 (22.6%)", y="PC1 (37.9%)") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="none"
  )


#ggsave(p, filename = "expPcaComp6.pdf", height = 9, width = 5)
########################

data("mDB")
output = TF_Selection(GSDB = mDB, DErslt = macrophageData$DErslt, minSize = 8, nperm = 1000, qval = .05, compList = compare_list, ntop = NULL, nameFile = "GSEArslts")

########### netact #################
tfs <- Reselect_TFs(GSEArslt = output$GSEArslt, qval = c(.05, .01, .01, .01, .01, .01), combine_TFs = T)
tfs_sep <- Reselect_TFs(GSEArslt = output$GSEArslt, qval = c(.05, .01, .01, .01, .01, .01), combine_TFs = F)
tfs_1 <- tfs_sep$`UT-IL4_2h`
tfs_2 <- tfs_sep$`UT-IFNg_2h`
tfs_3 <- tfs_sep$`UT-IFNg_IL4_2h`
tfs_4 <- tfs_sep$`UT-IL4_4h`
tfs_5 <- tfs_sep$`UT-IFNg_4h`
tfs_6 <- tfs_sep$`UT-IFNg_IL4_4h`

act.me <- TF_Activity(tfs, mDB, macrophageData$neweset, macrophageData$DErslt$Overall, with_weight = TRUE, useDatabaseSign = F, useCorSign = T, if_module = F)
acts_mat <- act.me$all_activities
# Stat6 activity sign was changed as we know it is wrong (using information from the paper)


acts_mat["Stat6",] <- - acts_mat["Stat6",]

####look up rational for flipping stat6

######### heatmap for initial - useless ######
exprs_mat = exprs(macrophageData$neweset)[rownames(acts_mat), ]
new_activity  <- acts_mat[,sort(colnames(acts_mat))]
data <- exprs_mat

n = 7

library(metafolio)
groupColors = gg_color_hue(n)

groupColors <- groupColors[c(7,5,1,3,6,2,4)]
# plot(1:n, pch = 16, cex = 2, col = groupColors)
#groupColors <- as.list(groupColors)
names(groupColors) <- plotData$type[1:7]

H1 =  ComplexHeatmap::Heatmap(row_norm(new_activity), col = circlize::colorRamp2(c(-2, 0, 2), c("green3", "white", "red")),
                              cluster_columns = F, cluster_rows = T, show_row_dend = FALSE, name = "Activity", column_title = "Activity",
                              row_names_gp = grid::gpar(fontsize = 15), column_names_gp = grid::gpar(fontsize = 15), clustering_method_rows = "ward.D2",
                              clustering_distance_rows = function(x) as.dist((1-cor(t(x), method = "spear"))/2)) 
# ,
#                              top_annotation = HeatmapAnnotation(data.frame(Treatment = substr(colnames(new_activity),1,(nchar(colnames(new_activity))-3))), col = list(Treatment = groupColors)))
gs = rownames(new_activity)
gc = colnames(new_activity)
H2 = ComplexHeatmap::Heatmap(row_norm(data[gs, gc]), col = circlize::colorRamp2(c(-2, 0, 2), c("green3", "white", "red")),
                             cluster_columns = F, cluster_rows = T, show_row_dend = FALSE,name = "Expression", column_title = "Expression",
                             row_names_gp = grid::gpar(fontsize = 15), column_names_gp = grid::gpar(fontsize = 15), clustering_method_rows = "ward.D2",
                             clustering_distance_rows = function(x) as.dist((1-cor(t(x), method = "spear"))/2))
#,
#                             top_annotation = HeatmapAnnotation(data.frame(Treatment = substr(colnames(new_activity),1,(nchar(colnames(new_activity))-3))), col = list(Treatment = groupColors)))
H = H1 + H2
draw(H, show_annotation_legend = FALSE)

########### NETWORK CONSTRUCTION ##########################

# Final Network - useCor = T, miTh = 0.5, DPI=F
miTh = 0.4
minMiTh = 0.5#.7

tfsSelected <- rownames(acts_mat)[which(rownames(acts_mat) %in% sort(unique(as.character(c(tfs_1,tfs_4)))))]
actsMat <- acts_mat[tfsSelected,c("UT_R1", "UT_R2", "UT_R3",  "IL4_2h_R1",  "IL4_2h_R2",  "IL4_2h_R3", "IL4_4h_R1", "IL4_4h_R2", "IL4_4h_R3"   )]

#maxInteractions = 500
tfLinks1 = TF_Filter(actsMat, mDB, miTh = miTh, nbins =3, maxTf = 150, 
                     maxInteractions = 500,  corMethod = "s", 
                     useCor = T, removeSignalling = FALSE, DPI = F, 
                     miDiff = 0.0, minMiTh = minMiTh)
dim(tfLinks1)
tfs_network <- sort(union(tfLinks1[,1],tfLinks1[,2]))
tfs_network
plot_network(tfLinks1)

tfsSelected <- rownames(acts_mat)[which(rownames(acts_mat) %in% sort(unique(as.character(c(tfs_2,tfs_5)))))]
actsMat <- acts_mat[tfsSelected,c("UT_R1", "UT_R2", "UT_R3",  "IFNg_2h_R1",  "IFNg_2h_R2",  "IFNg_2h_R3", "IFNg_4h_R1", "IFNg_4h_R2", "IFNg_4h_R3"   )]
tfLinks2 = TF_Filter(actsMat, mDB, miTh = miTh, nbins =3, maxTf = 150, 
                     maxInteractions = 500,  corMethod = "s", 
                     useCor = T, removeSignalling = FALSE, DPI = F, 
                     miDiff = 0.0, minMiTh = minMiTh)
dim(tfLinks2)
#plot_network(tfLinks2)

tfsSelected <- rownames(acts_mat)[which(rownames(acts_mat) %in% sort(unique(as.character(c(tfs_3,tfs_6)))))]
actsMat <- acts_mat[tfsSelected,c("UT_R1", "UT_R2", "UT_R3",  "IFNg_IL4_2h_R1",  "IFNg_IL4_2h_R2",  "IFNg_IL4_2h_R3", "IFNg_IL4_4h_R1", "IFNg_IL4_4h_R2", "IFNg_IL4_4h_R3"   )]
tfLinks3 = TF_Filter(actsMat, mDB, miTh = miTh, nbins =3, maxTf = 150, 
                     maxInteractions = 500,  corMethod = "s", 
                     useCor = T, removeSignalling = FALSE, DPI = F, 
                     miDiff = 0.0, minMiTh = minMiTh)
dim(tfLinks3)
# plot_network(tfLinks3)



tfLinks <- rbind(tfLinks1,tfLinks2, tfLinks3)
# plot_network(tfLinks)

tfs_network <- sort(union(tfLinks[,1],tfLinks[,2]))
tfs_network
nonDupTfLinks <- rep(TRUE, length = length(tfLinks[,1]))
TfLinksType <- tfLinks[,3]
# Fix the signs if there is discrepancy in sign
for(i in seq_along(tfLinks[,1])){
  tmp <- which(tfLinks[,1] == tfLinks[i,1] & tfLinks[,2] == tfLinks[i,2] )
  intCount <- c(0,0)
  if(length(tmp)>1){
    for(j in 2:length(tmp)){
      intCount[tfLinks[tmp[j],3]] <- intCount[tfLinks[tmp[j],3]] +1
      nonDupTfLinks[tmp[j]] <- FALSE
    }
    if(intCount[1] == intCount[2]){
      TfLinksType[tmp[1]] <- 2
    }
    else {
      TfLinksType[tmp[1]] <- which.max(intCount)
    }
  }
}
tfLinksN <- tfLinks[nonDupTfLinks,]
tfLinksN[,3] <- TfLinksType[nonDupTfLinks]


tfs_network <- sort(union(tfLinksN[,1],tfLinksN[,2]))
tfs_network

plot_network(tfLinksN)

#######simulate network with 2000 models ########

#tmp <- sRACIPE::sracipeSimulate(tfLinksN, plots = F, plotToFile = F, numModels = 2000)
# tmp <- sRACIPE::simulateGRC(tfLinksN, plots = TRUE, plotToFile = TRUE, nClusters = 7, numModels = 2000)
# tmp <- sRACIPE::plotData(tmp, plotToFile = TRUE, nClusters = 7)
#tmp2 <- sRACIPE::sracipeSimulate(tfLinksN, plots = T, plotToFile = F, nClusters = 7, numModels = 2000)
library(sRACIPE)
#tmp <- sRACIPE::sracipeSimulate(tfLinksN, plots = F, plotToFile = F, numModels = 10000)
#saveRDS(tmp, file = "mac.net.BRC")
#saveRDS(tmp, file = "mac.net.BRC.05")
#saveRDS(tmp, file = "mac.net.BRC.01.05")
#tmp <- readRDS("mac.net.BRC")
#tmp <- readRDS("mac.net.BRC.05")
#tmp <- readRDS("mac.net.BRC.01.05")
#gex <- log2(t(assay(tmp)))
#pca <- prcomp(gex, center = T, scale = T)

#ggplot(as.data.frame(pca$x), aes(x = PC1, y  =PC2)) + geom_point()



############ switch signs ###############

nodes <- unique(c(tfLinksN$Source, tfLinksN$Target))

#idenify DEs and non DEs
non.DEGs <- row.names(acts_mat)[row.names(acts_mat) %in% row.names(macrophageData$DErslt$Overall$table)[macrophageData$DErslt$Overall$table$adj.P.Val >= .05]]
non.in.net <- non.DEGs[non.DEGs %in% nodes]
DEGs <- row.names(acts_mat)[row.names(acts_mat) %in% row.names(macrophageData$DErslt$Overall$table)[macrophageData$DErslt$Overall$table$adj.P.Val < .05]]
de.in.net <- DEGs[DEGs %in% nodes]

#identify DE and non DE targets
de.in.net.targets <- vector(mode = "list", length = length(de.in.net))
for (i in 1:length(de.in.net)){
  de.in.net.targets[[i]] <- row.names(act.me$all_list[de.in.net[[i]]][[1]])
}
names(de.in.net.targets) <- de.in.net
#View(de.in.net.targets)

non.in.net.targets <- vector(mode = "list", length = length(non.in.net))
for (i in 1:length(non.in.net)){
  non.in.net.targets[[i]] <- row.names(act.me$all_list[non.in.net[[i]]][[1]])
}
names(non.in.net.targets) <- non.in.net
#View(non.in.net.targets)

#### now need to find overlap between targets

#make pairwise combos of two lists
a <- length(non.in.net)
b <- length(de.in.net)
C <- expand.grid(1:a, 1:b)  
combos.vec <- vector(length = nrow(C))
o.mat <- as.data.frame(matrix(nrow = length(combos.vec), ncol = 10))
colnames(o.mat) <- c("pair", "tf1.targets", "tf2.targets", "tf1.non", "tf2.non", "A","B","C","D", "P.value")

for (i in 1:nrow(C)){
  combos.vec[[i]] <-paste(C[i,][[1]], C[i,][[2]], sep = " ")
}

P.mat <- matrix(nrow = length(non.in.net), ncol = length(de.in.net))
rownames(P.mat) <- non.in.net
colnames(P.mat) <- de.in.net
ex.genes <- macrophageData$DErslt$Overall$table$genes

for (i in 1:length(combos.vec)){
  tf1 <- non.in.net[[as.integer(strsplit(combos.vec[[i]], split = " ")[[1]][1])]]
  tf2 <- de.in.net[[as.integer(strsplit(combos.vec[[i]], split = " ")[[1]][2])]]
  o.mat[i,1] <- paste0(tf1, " and ", tf2)
  tf1.targets <- non.in.net.targets[[tf1]]
  tf2.targets <- de.in.net.targets[[tf2]]
  tf1.non <- setdiff(ex.genes, tf1.targets)
  tf2.non <- setdiff(ex.genes, tf2.targets)
  A <- length(intersect(tf1.targets, tf2.targets))
  B <- length(intersect(tf1.targets, tf2.non))
  C <- length(intersect(tf1.non, tf2.targets))
  D <- length(intersect(tf1.non, tf2.non))
  fish.mat <- matrix(c(A,C,B,D), nrow = 2, ncol=2)
  #fish <- fisher.test(fish.mat)
  fish <- fisher.test(fish.mat, alternative = "greater")
  #fish <- fisher.test(fish.mat, alternative = "less")
  o.mat[i,2] <- length(tf1.targets)
  o.mat[i,3] <- length(tf2.targets)
  o.mat[i,4] <- length(tf1.non)
  o.mat[i,5] <- length(tf2.non)
  o.mat[i,6] <- A
  o.mat[i,7] <- B
  o.mat[i,8] <- C
  o.mat[i,9] <- D
  o.mat[i,10] <- fish$p.value
  P.mat[tf1,tf2] <- fish$p.value
}
#View(P.mat)
#im(P.mat)
# for each non DE find top 5 most signifcant P values
# first make a matrix where each row is a non DE and each column is a DE it pairs with and the value is the p value

P.mat <- t(P.mat)
P.mat <- as.data.frame(P.mat)

fish.list <- vector(mode = "list")
for(i in 1:length(colnames(P.mat))){
  tmp <-P.mat[,i]
  names(tmp) <-rownames(P.mat)
  tmp <- tmp[tmp < .05]
  fish.list[[i]] <-names(tmp[order(tmp)])
}
length(fish.list)
names(fish.list) <- colnames(P.mat)
#View(fish.list)


#fish.list$Clock

#switch signs that DONT agree 
NDE <- "Clock"
DE <-"Arntl"
sign(acts_mat[NDE,])
sign(acts_mat[DE,])
if(sum(sign(acts_mat[NDE,])!= sign(acts_mat[DE,])) > (length(sign(acts_mat[NDE,])!= sign(acts_mat[DE,]))/2)){
  acts_mat[NDE,] <- -acts_mat[NDE,]
  print("switched")
}

######### qval .05 ##########


#switch signs that DONT agree 
NDE <- "Crebbp"
DE <-"Hif1a"
sign(acts_mat[NDE,])
sign(acts_mat[DE,])
if(sum(sign(acts_mat[NDE,]) != sign(acts_mat[DE,])) > (length(sign(acts_mat[NDE,])!= sign(acts_mat[DE,]))/2)){
  acts_mat[NDE,] <- -acts_mat[NDE,]
  print("switched")
}




########### NETWORK CONSTRUCTION 2 ##########################

# Final Network - useCor = T, miTh = 0.5, DPI=F
miTh = 0.4
minMiTh = 0.5#.7



tfsSelected <- rownames(acts_mat)[which(rownames(acts_mat) %in% sort(unique(as.character(c(tfs_1,tfs_4)))))]
actsMat <- acts_mat[tfsSelected,c("UT_R1", "UT_R2", "UT_R3",  "IL4_2h_R1",  "IL4_2h_R2",  "IL4_2h_R3", "IL4_4h_R1", "IL4_4h_R2", "IL4_4h_R3"   )]
tfLinks1 = TF_Filter(actsMat, mDB, miTh = miTh, nbins =3, maxTf = 150, 
                     maxInteractions = 500,  corMethod = "s", 
                     useCor = T, removeSignalling = FALSE, DPI = F, 
                     miDiff = 0.0, minMiTh = minMiTh)
dim(tfLinks1)
tfs_network <- sort(union(tfLinks1[,1],tfLinks1[,2]))
tfs_network
#plot_network(tfLinks1)

tfsSelected <- rownames(acts_mat)[which(rownames(acts_mat) %in% sort(unique(as.character(c(tfs_2,tfs_5)))))]
actsMat <- acts_mat[tfsSelected,c("UT_R1", "UT_R2", "UT_R3",  "IFNg_2h_R1",  "IFNg_2h_R2",  "IFNg_2h_R3", "IFNg_4h_R1", "IFNg_4h_R2", "IFNg_4h_R3"   )]
tfLinks2 = TF_Filter(actsMat, mDB, miTh = miTh, nbins =3, maxTf = 150, 
                     maxInteractions = 500,  corMethod = "s", 
                     useCor = T, removeSignalling = FALSE, DPI = F, 
                     miDiff = 0.0, minMiTh = minMiTh)
dim(tfLinks2)
#plot_network(tfLinks2)

tfsSelected <- rownames(acts_mat)[which(rownames(acts_mat) %in% sort(unique(as.character(c(tfs_3,tfs_6)))))]
actsMat <- acts_mat[tfsSelected,c("UT_R1", "UT_R2", "UT_R3",  "IFNg_IL4_2h_R1",  "IFNg_IL4_2h_R2",  "IFNg_IL4_2h_R3", "IFNg_IL4_4h_R1", "IFNg_IL4_4h_R2", "IFNg_IL4_4h_R3"   )]
tfLinks3 = TF_Filter(actsMat, mDB, miTh = miTh, nbins =3, maxTf = 150, 
                     maxInteractions = 500,  corMethod = "s", 
                     useCor = T, removeSignalling = FALSE, DPI = F, 
                     miDiff = 0.0, minMiTh = minMiTh)
dim(tfLinks3)
# plot_network(tfLinks3)



tfLinks <- rbind(tfLinks1,tfLinks2, tfLinks3)
# plot_network(tfLinks)

tfs_network <- sort(union(tfLinks[,1],tfLinks[,2]))
tfs_network
nonDupTfLinks <- rep(TRUE, length = length(tfLinks[,1]))
TfLinksType <- tfLinks[,3]
# Fix the signs if there is discrepancy in sign
for(i in seq_along(tfLinks[,1])){
  tmp <- which(tfLinks[,1] == tfLinks[i,1] & tfLinks[,2] == tfLinks[i,2] )
  intCount <- c(0,0)
  if(length(tmp)>1){
    for(j in 2:length(tmp)){
      intCount[tfLinks[tmp[j],3]] <- intCount[tfLinks[tmp[j],3]] +1
      nonDupTfLinks[tmp[j]] <- FALSE
    }
    if(intCount[1] == intCount[2]){
      TfLinksType[tmp[1]] <- 2
    }
    else {
      TfLinksType[tmp[1]] <- which.max(intCount)
    }
  }
}
tfLinksN <- tfLinks[nonDupTfLinks,]
tfLinksN[,3] <- TfLinksType[nonDupTfLinks]





######### for qval .05 only ##########

tfLinksN <- tfLinksN[-(which(tfLinksN$Source == "Bmi1" | tfLinksN$Target == "Bmi1")),]
tfLinksN <- tfLinksN[-(which(tfLinksN$Source == "Hdac3" | tfLinksN$Target == "Hdac3")),]
tfLinksN <- tfLinksN[-(which(tfLinksN$Source == "Tfdp1" | tfLinksN$Target == "Tfdp1")),]





plot_network(tfLinksN)
tfs_network <- sort(union(tfLinksN[,1],tfLinksN[,2]))


############ heatmaps again ################
new_activity  <- acts_mat[,sort(colnames(acts_mat))]
data <- exprs_mat

n = 7

library(metafolio)
groupColors = gg_color_hue(n)

groupColors <- groupColors[c(7,5,1,3,6,2,4)]
# plot(1:n, pch = 16, cex = 2, col = groupColors)
#groupColors <- as.list(groupColors)
names(groupColors) <- plotData$type[1:7]

H1 =  ComplexHeatmap::Heatmap(row_norm(new_activity), col = circlize::colorRamp2(c(-2, 0, 2), c("blue3", "white", "red")),
                              cluster_columns = F, cluster_rows = T, show_row_dend = FALSE, name = "Activity", column_title = "Activity",
                              row_names_gp = grid::gpar(fontsize = 10), column_names_gp = grid::gpar(fontsize = 7), clustering_method_rows = "ward.D2",
                              clustering_distance_rows = function(x) as.dist((1-cor(t(x), method = "spear"))/2)) 
# ,
#                              top_annotation = HeatmapAnnotation(data.frame(Treatment = substr(colnames(new_activity),1,(nchar(colnames(new_activity))-3))), col = list(Treatment = groupColors)))
gs = rownames(new_activity)
gc = colnames(new_activity)
H2 = ComplexHeatmap::Heatmap(row_norm(data[gs, gc]), col = circlize::colorRamp2(c(-2, 0, 2), c("blue3", "white", "red")),
                             cluster_columns = F, cluster_rows = T, show_row_dend = FALSE,name = "Expression", column_title = "Expression",
                             row_names_gp = grid::gpar(fontsize = 10), column_names_gp = grid::gpar(fontsize = 7), clustering_method_rows = "ward.D2",
                             clustering_distance_rows = function(x) as.dist((1-cor(t(x), method = "spear"))/2))
#,
#                             top_annotation = HeatmapAnnotation(data.frame(Treatment = substr(colnames(new_activity),1,(nchar(colnames(new_activity))-3))), col = list(Treatment = groupColors)))
H = H1 + H2
draw(H, show_annotation_legend = FALSE)



############ simulate new network with only 2k genes ##########
#length(unique(tfLinksN$Source,tfLinksN$Target))
#new.net.rset <- sRACIPE::sracipeSimulate(tfLinksN, plots = F, plotToFile = F, numModels = 500)
#save(new.net.rset, file = "new.net.rset.05.01_corrected")
#load("new.net.rset")
load("new.net.rset.05.01_corrected")


gex <- log2(t(assay(new.net.rset)))
pca <- prcomp(gex, center = T, scale = T)

ggplot(as.data.frame(pca$x), aes(x = PC1, y  =PC2)) + geom_point()

#t <- getwd()
#setwd("/Users/c-clausb/Downloads")
#sracipePlotData(new.net.rset, nClusters = 7, umapPlot = T)
#setwd(t)

k <- kmeans(pca$x, centers = 7)
#save(k, file = "/Users/c-clausb/Desktop/Netact/reneedsomehelpfornetactmodeling/macroK_corrected")
#load("/Users/c-clausb/Desktop/Netact/reneedsomehelpfornetactmodeling/macroK_corrected")
ggplot(as.data.frame(pca$x), aes(x = PC1, y  =PC2, color = as.factor(k$cluster))) + geom_point()

#"Creb3" %in% unique(c(tfLinksN$Source,tfLinksN$Target))
#"Atf2" %in% unique(c(tfLinksN$Source,tfLinksN$Target))
#"Ncoa1" %in% unique(c(tfLinksN$Source,tfLinksN$Target))
#"Creb5" %in% unique(c(tfLinksN$Source,tfLinksN$Target))
#"Atf6b" %in% unique(c(tfLinksN$Source,tfLinksN$Target))
#"Dach1" %in% unique(c(tfLinksN$Source,tfLinksN$Target))
#"Pelp1" %in% unique(c(tfLinksN$Source,tfLinksN$Target))
#"Smad1" %in% unique(c(tfLinksN$Source,tfLinksN$Target))
#"Myc" %in% unique(c(tfLinksN$Source,tfLinksN$Target))


############ mapping back to experimetnal data  using expression data###########

#View(combat_edata)
#colnames(combat_edata)
combat_edata.s <- t(scale(t(combat_edata), center = T, scale = T))
#View(combat_edata.s)
group.1 <- data.frame( MeanGex = rowMeans(combat_edata.s[,colnames(combat_edata.s)[grep(pattern = "UT_R", x = colnames(combat_edata.s))]]))
group.2 <- data.frame(MeanGex = rowMeans(combat_edata.s[,colnames(combat_edata.s)[grep(pattern = "IFNg_2h", x = colnames(combat_edata.s))]]))
group.3 <- data.frame(MeanGex = rowMeans(combat_edata.s[,colnames(combat_edata.s)[grep(pattern = "IFNg_4h", x = colnames(combat_edata.s))]]))
group.4 <- data.frame(MeanGex = rowMeans(combat_edata.s[,colnames(combat_edata.s)[grep(pattern = "IFNg_IL4_2h", x = colnames(combat_edata.s))]]))
group.5 <- data.frame(MeanGex = rowMeans(combat_edata.s[,colnames(combat_edata.s)[grep(pattern = "IFNg_IL4_4h", x = colnames(combat_edata.s))]]))
group.6 <- data.frame(MeanGex = rowMeans(combat_edata.s[,colnames(combat_edata.s)[grep(pattern = "IL4_2h", x = colnames(combat_edata.s))]]))
group.7 <- data.frame(MeanGex = rowMeans(combat_edata.s[,colnames(combat_edata.s)[grep(pattern = "IL4_4h", x = colnames(combat_edata.s))]]))

ex.genes <- rownames(combat_edata.s)

map.me <- t(scale(log2(t(assay(new.net.rset))),center = T, scale = T))
#class(map.me)
#dim(map.me)
#rownames(map.me)
sim.genes <- rownames(map.me)

sub.1 <- group.1[sim.genes,]
names(sub.1) <- sim.genes
sub.2 <- group.2[sim.genes,]
names(sub.2) <- sim.genes
sub.3 <- group.3[sim.genes,]
names(sub.3) <- sim.genes
sub.4 <- group.4[sim.genes,]
names(sub.4) <- sim.genes
sub.5 <- group.5[sim.genes,]
names(sub.5) <- sim.genes
sub.6 <- group.6[sim.genes,]
names(sub.6) <- sim.genes
sub.7 <- group.7[sim.genes,]
names(sub.7) <- sim.genes


#all.sub <- t(combat_edata.s[sim.genes,])
#all.sub.pca <- prcomp(all.sub, center = T)
#rownames(all.sub.pca$x)


rownames(all.sub.pca$x)[grep(pattern = "UT_", rownames(all.sub.pca$x))] <- "UT"
rownames(all.sub.pca$x)[grep(pattern = "IFNg_IL4_2h_R", rownames(all.sub.pca$x))] <- "IFNg_IL4_2h"
rownames(all.sub.pca$x)[grep(pattern = "IFNg_IL4_4h_R", rownames(all.sub.pca$x))] <- "IFNg_IL4_4h"
rownames(all.sub.pca$x)[grep(pattern = "IL4_2h_R", rownames(all.sub.pca$x))] <- "IL4_2h"
rownames(all.sub.pca$x)[grep(pattern = "IL4_4h_R", rownames(all.sub.pca$x))] <- "IL4_4h"
rownames(all.sub.pca$x)[grep(pattern = "IFNg_2h_R", rownames(all.sub.pca$x))] <- "IFNg_2h"
rownames(all.sub.pca$x)[grep(pattern = "IFNg_4h_R", rownames(all.sub.pca$x))] <- "IFNg_4h"


#ggplot(as.data.frame(all.sub.pca$x), aes (x = PC1, y = -PC2, color = as.factor(rownames(all.sub.pca$x)))) + geom_point()

#View(map.me)
pca.tmp <- prcomp(t(map.me), center = F, scale  = F)
ggplot(as.data.frame(pca.tmp$x), aes(x = PC1, y = PC2)) + geom_point()

dmat <- matrix(nrow = 500, ncol = 7)
colnames(dmat) <- c("UT", "IFNg_2h", "IFNg_4h", "IFNg_IL4_2h", "IFNg_IL4_4h", "IL4_2h", "IL4_4h")

euclidean <- function(a, b) sqrt(sum((a - b)^2))
#euclidean <- function(a, b) sqrt((a - b) %*% (a - b))

for(i in 1:ncol(map.me)){
  dmat[i, 1] <- euclidean(map.me[,i], sub.1 )
  dmat[i, 2] <- euclidean(map.me[,i], sub.2 )
  dmat[i, 3] <- euclidean(map.me[,i], sub.3 )
  dmat[i, 4] <- euclidean(map.me[,i], sub.4 )
  dmat[i, 5] <- euclidean(map.me[,i], sub.5 )
  dmat[i, 6] <- euclidean(map.me[,i], sub.6 )
  dmat[i, 7] <- euclidean(map.me[,i], sub.7 )
}

View(dmat)

model.assignment <- integer(length = 500)
for(i in 1:nrow(dmat)){
  model.assignment[[i]] <- colnames(dmat)[which.min(dmat[i,])]
}

View(model.assignment)
ggplot(as.data.frame(pca.tmp$x), aes(x = PC1, y = PC2, color = model.assignment)) + geom_point()
ggplot(as.data.frame(pca.tmp$x), aes(x = PC1, y = PC2, color = model.assignment)) + geom_point() +stat_ellipse()

############ mapping back to experimental data  using activity data###########




act.mat <- t(act.me$all_activities)
act.mat.s <- t(scale(act.mat))



#act.mat.s

#colnames(act.mat.s)

group.1 <- data.frame( MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "UT_R", x = colnames(act.mat.s))]]))
group.2 <- data.frame(MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "IFNg_2h", x = colnames(act.mat.s))]]))
group.3 <- data.frame(MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "IFNg_4h", x = colnames(act.mat.s))]]))
group.4 <- data.frame(MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "IFNg_IL4_2h", x = colnames(act.mat.s))]]))
group.5 <- data.frame(MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "IFNg_IL4_4h", x = colnames(act.mat.s))]]))
group.6 <- data.frame(MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "IL4_2h", x = colnames(act.mat.s))]]))
group.7 <- data.frame(MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "IL4_4h", x = colnames(act.mat.s))]]))




sub.1 <- group.1[sim.genes,]
names(sub.1) <- sim.genes
sub.2 <- group.2[sim.genes,]
names(sub.2) <- sim.genes
sub.3 <- group.3[sim.genes,]
names(sub.3) <- sim.genes
sub.4 <- group.4[sim.genes,]
names(sub.4) <- sim.genes
sub.5 <- group.5[sim.genes,]
names(sub.5) <- sim.genes
sub.6 <- group.6[sim.genes,]
names(sub.6) <- sim.genes
sub.7 <- group.7[sim.genes,]
names(sub.7) <- sim.genes

all.sub <- act.mat.s[sim.genes,]
all.sub.2 <- act.mat[,sim.genes]
all.sub.pca <- prcomp(t(all.sub), center = F)
#rownames(all.sub.pca$x)


rownames(all.sub.pca$x)[grep(pattern = "UT_", rownames(all.sub.pca$x))] <- "UT"
rownames(all.sub.pca$x)[grep(pattern = "IFNg_IL4_2h_R", rownames(all.sub.pca$x))] <- "IFNg_IL4_2h"
rownames(all.sub.pca$x)[grep(pattern = "IFNg_IL4_4h_R", rownames(all.sub.pca$x))] <- "IFNg_IL4_4h"
rownames(all.sub.pca$x)[grep(pattern = "IL4_2h_R", rownames(all.sub.pca$x))] <- "IL4_2h"
rownames(all.sub.pca$x)[grep(pattern = "IL4_4h_R", rownames(all.sub.pca$x))] <- "IL4_4h"
rownames(all.sub.pca$x)[grep(pattern = "IFNg_2h_R", rownames(all.sub.pca$x))] <- "IFNg_2h"
rownames(all.sub.pca$x)[grep(pattern = "IFNg_4h_R", rownames(all.sub.pca$x))] <- "IFNg_4h"

########## plot experimental activity of network genes #########
summary(all.sub.pca)

p <- ggplot(as.data.frame(all.sub.pca$x), aes(x=-PC2, y=PC1, color = as.factor(rownames(all.sub.pca$x)))) +
  geom_point(size=5) + 
  labs(y="PC1 (81.2%)", x="PC2 (11.8%)", col="Treatment", shape="Batch") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="none"
    )
p


plot.me.2 <- all.sub.pca$x[,1:2]

plot.me.2 <- cbind(plot.me.2, rownames(all.sub.pca$x), c(0,2,2,2,4,4,4,0,2,2,2,4,4,4,0,2,2,2,4,4,4))
colnames(plot.me.2)<- c("PC1", "PC2", "Treatment", "Hour")
plot.me.2 <- as.data.frame(plot.me.2)
plot.me.2$PC1 <- as.numeric(as.character(plot.me.2$PC1))
plot.me.2$PC2 <- as.numeric(as.character(plot.me.2$PC2))
plot.me.2$Treatment <- as.character(plot.me.2$Treatment)


plot.me.2$Treatment[plot.me.2$Treatment == "IL4_2h" |plot.me.2$Treatment == "IL4_4h"] <- "IL4"
plot.me.2$Treatment[plot.me.2$Treatment == "IFNg_2h" |plot.me.2$Treatment == "IFNg_4h"] <- "IFNg"
plot.me.2$Treatment[plot.me.2$Treatment == "IFNg_IL4_2h" |plot.me.2$Treatment == "IFNg_IL4_4h"] <- "IFNg_IL4"



plot.me.2$Treatment <- as.factor(plot.me.2$Treatment)
plot.me.2$Hour <- as.factor(plot.me.2$Hour)

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")                      
ggplot(plot.me.2, aes(x=-PC2, y=PC1, col=Treatment, shape=Hour)) +
  geom_point(size=5) + 
  labs(y="PC1 (81.2%)", x="PC2 (11.8%)", col="Treatment", shape="Hour") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    legend.position="none")+
  scale_colour_manual(values=cbbPalette)







ggplot(as.data.frame(all.sub.pca$x), aes (x = PC1, y = -PC2, color = as.factor(rownames(all.sub.pca$x)))) + geom_point()

#View(map.me)

ggplot(as.data.frame(pca.tmp$x), aes(x = PC1, y = PC2)) + geom_point()






dmat <- matrix(nrow = 500, ncol = 7)
colnames(dmat) <- c("UT", "IFNg_2h", "IFNg_4h", "IFNg_IL4_2h", "IFNg_IL4_4h", "IL4_2h", "IL4_4h")

euclidean <- function(a, b) sqrt(sum((a - b)^2))
#euclidean <- function(a, b) sqrt((a - b) %*% (a - b))

for(i in 1:ncol(map.me)){
  dmat[i, 1] <- euclidean(map.me[,i], sub.1 )
  dmat[i, 2] <- euclidean(map.me[,i], sub.2 )
  dmat[i, 3] <- euclidean(map.me[,i], sub.3 )
  dmat[i, 4] <- euclidean(map.me[,i], sub.4 )
  dmat[i, 5] <- euclidean(map.me[,i], sub.5 )
  dmat[i, 6] <- euclidean(map.me[,i], sub.6 )
  dmat[i, 7] <- euclidean(map.me[,i], sub.7 )
}

#View(dmat)

model.assignment <- integer(length = 500)
for(i in 1:nrow(dmat)){
  model.assignment[[i]] <- colnames(dmat)[which.min(dmat[i,])]
}

View(model.assignment)

ggplot(as.data.frame(pca.tmp$x), aes(x = PC1, y = PC2, color = model.assignment)) + geom_point()

ggplot(as.data.frame(pca.tmp$x), aes(x = PC1, y = PC2, color = model.assignment)) + geom_point() + stat_ellipse()



#View(model.assignment)
model.assignment[model.assignment == "IFNg_IL4_2h" | model.assignment == "IFNg_IL4_4h"] <- "IFNg_IL4"
model.assignment[model.assignment == "IFNg_2h" | model.assignment == "IFNg_4h"] <- "IFNg"
model.assignment[model.assignment == "IL4_2h" | model.assignment == "IL4_4h"] <- "IL4"

ggplot(as.data.frame(pca.tmp$x), aes(x = PC1, y = PC2, color = model.assignment)) + geom_point()
ggplot(as.data.frame(pca.tmp$x), aes(x = PC1, y = PC2, color = model.assignment)) + geom_point() + stat_ellipse()

circuit <- sracipeCircuit(new.net.rset)
#save(circuit, file = "macro_circ_corrected")

########## final simulation #######


load("/Users/c-clausb/Downloads/mac.final.gex.corrected")
load("/Users/c-clausb/Downloads/mac.final.circ.corrected")

act.mat <- t(act.me$all_activities)
act.mat.s <- t(scale(act.mat))
sim.genes <- rownames(gex)


#act.mat.s

#colnames(act.mat.s)

group.1 <- data.frame( MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "UT_R", x = colnames(act.mat.s))]]))
group.2 <- data.frame(MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "IFNg_2h", x = colnames(act.mat.s))]]))
group.3 <- data.frame(MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "IFNg_4h", x = colnames(act.mat.s))]]))
group.4 <- data.frame(MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "IFNg_IL4_2h", x = colnames(act.mat.s))]]))
group.5 <- data.frame(MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "IFNg_IL4_4h", x = colnames(act.mat.s))]]))
group.6 <- data.frame(MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "IL4_2h", x = colnames(act.mat.s))]]))
group.7 <- data.frame(MeanGex = rowMeans(act.mat.s[,colnames(act.mat.s)[grep(pattern = "IL4_4h", x = colnames(act.mat.s))]]))




sub.1 <- group.1[sim.genes,]
names(sub.1) <- sim.genes
sub.2 <- group.2[sim.genes,]
names(sub.2) <- sim.genes
sub.3 <- group.3[sim.genes,]
names(sub.3) <- sim.genes
sub.4 <- group.4[sim.genes,]
names(sub.4) <- sim.genes
sub.5 <- group.5[sim.genes,]
names(sub.5) <- sim.genes
sub.6 <- group.6[sim.genes,]
names(sub.6) <- sim.genes
sub.7 <- group.7[sim.genes,]
names(sub.7) <- sim.genes

load("/Users/c-clausb/Desktop/Netact/object/gex")
load("/Users/c-clausb/Desktop/Netact/object/circ")
#save(gex, file ="/Users/c-clausb/Desktop/Netact/object/gex" )
#save(circ, file = "/Users/c-clausb/Desktop/Netact/object/circ")


gex <- log2(t(gex))

pca <- prcomp(gex, center = T, scale =T)
pca.sim <- pca
ggplot(as.data.frame(pca$x), aes(x = -PC1, y = PC2)) + geom_point()
map.me <- t(scale(gex,center = T, scale = T))
dmat <- matrix(nrow = 10000, ncol = 7)
colnames(dmat) <- c("UT", "IFNg_2h", "IFNg_4h", "IFNg_IL4_2h", "IFNg_IL4_4h", "IL4_2h", "IL4_4h")

euclidean <- function(a, b) sqrt(sum((a - b)^2))
#euclidean <- function(a, b) sqrt((a - b) %*% (a - b))

for(i in 1:ncol(map.me)){
  dmat[i, 1] <- euclidean(map.me[,i], sub.1 )
  dmat[i, 2] <- euclidean(map.me[,i], sub.2 )
  dmat[i, 3] <- euclidean(map.me[,i], sub.3 )
  dmat[i, 4] <- euclidean(map.me[,i], sub.4 )
  dmat[i, 5] <- euclidean(map.me[,i], sub.5 )
  dmat[i, 6] <- euclidean(map.me[,i], sub.6 )
  dmat[i, 7] <- euclidean(map.me[,i], sub.7 )
}

#View(dmat)

model.assignment <- integer(length = 10000)
for(i in 1:nrow(dmat)){
  model.assignment[[i]] <- colnames(dmat)[which.min(dmat[i,])]
}

#View(model.assignment)

ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = model.assignment)) + geom_point()
summary(pca)

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = model.assignment)) +
  geom_point(size=5) + 
  labs(x="PC1 (36.6.6%)", y="PC2 (12.4%)", col="Treatment") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    legend.position="none"
  )+
  scale_colour_manual(values=cbbPalette)








#ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = model.assignment)) + geom_point() + stat_ellipse()



#View(model.assignment)
model.assignment[model.assignment == "IFNg_IL4_2h" | model.assignment == "IFNg_IL4_4h"] <- "IFNg_IL4"
model.assignment[model.assignment == "IFNg_2h" | model.assignment == "IFNg_4h"] <- "IFNg"
model.assignment[model.assignment == "IL4_2h" | model.assignment == "IL4_4h"] <- "IL4"

#save(model.assignment, file = "/Users/c-clausb/Downloads/MacFinalModelAssignment")
load("/Users/c-clausb/Downloads/MacFinalModelAssignment")


ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = model.assignment)) + geom_point()

summary(pca)




ggplot(as.data.frame(pca$x), aes(x = -PC1, y = PC2, color = model.assignment)) +
  geom_point(size=2) + 
  labs(x="PC1 (36.3%)", y="PC2 (12.3%)", col="Treatment") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    legend.position="none"
  ) 



#ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = model.assignment)) + geom_point() + stat_ellipse()

########## logistic regression to identify number of models in transition #######

load("/Users/c-clausb/Downloads/mac.final.gex.corrected")
load("/Users/c-clausb/Downloads/mac.final.circ.corrected")
gex <- log2(t(gex))
pca <- prcomp(gex, center = T, scale =T)
ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2)) + geom_point()
tmp.k <- kmeans(gex, centers = 30)
#save(tmp.k, file = "mac.tmp.k")
load("mac.tmp.k")
tmp.k <- tmp.k$cluster  
ggplot(as.data.frame(pca$x), aes(x = -PC1, y = PC2, color =as.factor(tmp.k) )) + geom_point() +
  theme_classic() +
  theme(legend.position="none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))

tmp <- tmp.k
tmp[tmp != 17] <- 2
ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color =as.factor(tmp) )) + geom_point()

i.1 <- which(tmp.k == 2 |tmp.k == 3 |tmp.k == 7 |tmp.k == 8 | tmp.k == 9| tmp.k == 10 | tmp.k == 12| tmp.k == 13| tmp.k == 14| tmp.k == 19|tmp.k == 21 |tmp.k == 22|tmp.k == 23|tmp.k == 24|tmp.k == 26|tmp.k == 27|tmp.k == 28|tmp.k == 30)
i.2 <- which(tmp.k == 1 |tmp.k == 4  |tmp.k == 29|tmp.k == 16|tmp.k == 5  )
i.3 <- which(tmp.k == 6 |tmp.k == 15 |tmp.k == 18| tmp.k == 20|tmp.k == 25|tmp.k == 11|tmp.k == 17 )
tmp.k.2 <- tmp.k
tmp.k.2[i.1] <- 1
tmp.k.2[i.2] <- 2
tmp.k.2[i.3] <- 3

summary(pca)

ggplot(as.data.frame(pca$x), aes(x = -PC1, y = PC2, color =as.factor(tmp.k.2) )) + geom_point() +
  theme_classic() +
  labs(x="PC1 (36.3%)", y="PC2 (12.3%)")+
  theme(legend.position="none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) + 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
################


lmpca <- lm(PC2~PC1, data =as.data.frame(pca$x[i.2,]) )
plot(as.data.frame(pca$x[i.2,1:2]))
abline(lmpca)
summary(lmpca)
coef(lmpca)
#save(lmpca, file = "mac_lmpca")
load("mac_lmpca")

tmp.pca <- as.data.frame(pca$x)
tmp.pca$PC1 <- -(tmp.pca$PC1)
lmpca <- lm(PC2~PC1, data =tmp.pca[i.2,] )


ggplot(as.data.frame(pca$x[i.2,]), aes(x = -PC1, y = PC2)) +
  geom_point() +
  theme_classic()

ggplot(as.data.frame(pca$x[i.2,]), aes(x = -PC1, y = PC2, color = model.assignment[i.2])) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))



ggplot(as.data.frame(pca$x[i.2,]), aes(x = -PC1, y = PC2)) +
  geom_point() +
  geom_abline(intercept = 0.2498239, slope = -0.8129740, col = "red") +
  theme_classic()



plot(as.data.frame(pca$x)[i.2,1:2])
abline(lmpca)

my.model <-function(point){
  return( (point*-0.8129740)+0.2498239)
}

model.line <- unlist(lapply(X = seq(from = -6, to = 3, by = 0.1), FUN = my.model))
line.mat <- cbind(seq(from = -6, to = 3, by = 0.1), model.line)
plot(as.data.frame(tmp.pca)[i.2,1:2])



plot(as.data.frame(tmp.pca)[i.2,1:2])
lines(line.mat[,1], line.mat[,2], col = 'red', lwd = 1)





nrow(line.mat)
91/3
121-82
#2 bins of 41 and one of 39

bin.1 <- line.mat[1:30,]
bin.2 <- line.mat[31:60,]
bin.3 <- line.mat[61:91,]
nrow(bin.3)


ggplot(as.data.frame(pca$x)[i.3,1:2], aes(x = PC1, y = PC2)) + 
  geom_point() +
  geom_segment(aes(x = bin.1[1,1], y = bin.1[1,2],xend =bin.1[30,1], yend =bin.1[30,2])) +
  geom_segment(aes(x = bin.2[1,1], y = bin.2[1,2],xend =bin.2[30,1], yend =bin.2[30,2])) +
  geom_segment(aes(x = bin.3[1,1], y = bin.3[1,2],xend =bin.3[31,1], yend =bin.3[31,2]))

plot(as.data.frame(pca$x)[i.3,1:2])
lines(bin.1[,1], bin.1[,2], col = 'green', lwd = 1)
lines(bin.2[,1], bin.2[,2], col = 'red', lwd = 1)
lines(bin.3[,1], bin.3[,2], col = 'blue', lwd = 1)


dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  return(d)
} 

pca$x[i.3[[1]],1:2]

dist.i.3 <- matrix(nrow = length(i.3), ncol = 3)

for(i in 1:nrow(dist.i.3)){
  a <- pca$x[i.3[[i]],1:2]
  dist.i.3[i,1] <- dist2d(a,bin.1[1,], bin.1[30,])
  dist.i.3[i,2] <- dist2d(a,bin.2[1,], bin.2[30,])
  dist.i.3[i,3] <- dist2d(a,bin.3[1,], bin.3[31,])
}

dist.i.3.final <- vector(length = length(i.3))

for(i in 1:length(dist.i.3.final)){
  dist.i.3.final[[i]] <- which.min(dist.i.3[i,])
}

dist.i.3.final

############ plot PCA + lines + assignment ##########


final.plot.mat <- as.data.frame(matrix(nrow = 10000, ncol = 3))
colnames(final.plot.mat) <- c("PC1", "PC2", "Assignment")
final.plot.mat$PC1 <- pca$x[,1]
final.plot.mat$PC2 <- pca$x[,2]
final.plot.mat$Assignment <- model.assignment

final.plot.mat[i.3,3] <- dist.i.3.final


ggplot(data = final.plot.mat[i.3,], aes(x = PC1, y = PC2, color = as.factor(Assignment))) + geom_point()

ggplot(data = final.plot.mat[i.3,], aes(x = PC1, y = PC2, color = as.factor(Assignment))) + geom_point()

hist(final.plot.mat[,1:2])



########### project points to linear regression line ############

library(LearnGeom)
ProjectPoint()



p1 <- tmp.pca[i.2,1:2][1,]
tmp.pca[i.2,1:2][1,]


summary(lmpca)
line.start <- c(-6,5.1276679)
line.end <-c(3,-3.87)

line <- CreateLinePoints(line.start,line.end)

ProjectPoint(p1, line)

projected <- as.data.frame(matrix(nrow = nrow(pca$x[i.3,1:2]), ncol = 3))

colnames(projected) <- c("X", "Y", "Assignment")

for(i in 1:nrow(pca$x[i.3,1:2])){
  p1 <- pca$x[i.3,1:2][i,]
  projected[i,1] <- ProjectPoint(p1, line)[[1]]
  projected[i,2] <- ProjectPoint(p1, line)[[2]]
  projected[i,3] <- model.assignment[i.3][i]
}

#View(projected)

ggplot(as.data.frame(pca$x[i.3,1:2]), aes(x = PC1, y = PC2))+geom_point()

ggplot(projected, aes(x = X, y = Y, color = model.assignment[i.3])) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))

projected.IFNg <- projected[projected$Assignment == "IFNg",]
projected.IL4 <- projected[projected$Assignment == "IL4",]
projected.IFNg_IL4 <- projected[projected$Assignment == "IFNg_IL4",]


hist(projected.IFNg[,1], col = "red")
hist(projected.IL4[,1], col = "blue")
hist(projected.IFNg_IL4[,1], col = "green")
hist(projected.IL4[,1], col = "blue", add = T)
hist(projected.IFNg_IL4[,1], col = "green", add = T)

projected.all <- rbind(projected.IFNg,projected.IL4,projected.IFNg_IL4 )

ggplot(projected, aes(x=X, fill=Assignment)) +
  geom_histogram( color='#e9ecef', alpha=0.6, position='identity', bins = 15) 

ggplot(projected, aes(x=X, fill=Assignment)) +
  geom_density( color='#e9ecef', alpha=0.6, position='identity') +
  xlab("PC1") +
  ylab("Density") + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))


























########## making arbitrary lines to identify transitory models ###########

mk.line <- function(m, x, b){
  y <- (m*x)+b
  return(y)
}

line.1.x <- seq(from = (2), to = 9, length.out = 10)
line.1.y <- unlist(lapply(line.1.x, mk.line, m = (12/6), b = -11 ))
line.2.x <- seq(from = (-3), to = 4, length.out = 10)
line.2.y <- unlist(lapply(line.2.x, mk.line, m = (12/6), b = 0 ))


line.1.df <- as.data.frame(cbind(line.1.x, line.1.y))
line.2.df <- as.data.frame(cbind(line.2.x, line.2.y))
colnames(line.1.df) <-c("X", "Y")

ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = model.assignment)) +
  geom_point(size=2) + 
  labs(x="PC1 (36.3%)", y="PC2 (12.3%)", col="Treatment") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    legend.position="none") +
  geom_segment(aes(x = line.1.df[1,1], y = line.1.df[1,2], xend = line.1.df[nrow(line.1.df),1], yend = line.1.df[nrow(line.1.df),2] )) +
  geom_segment(aes(x = line.2.df[1,1], y = line.2.df[1,2], xend = line.2.df[nrow(line.2.df),1], yend = line.2.df[nrow(line.2.df),2] )) +
  scale_colour_manual(values=cbbPalette)
 
find.y <- function(m, x, b){
  return((m*x) +b)
}
find.loc <- function(x,y){
  y.1 <- find.y(12/6, x, -11)
  y.2 <-find.y(12/6, x, 0)
  if (y>y.1 & y < y.2){
    return("mid")
  } else if (y > y.2){
    return("upper")
  } else if (y < y.1){
    return("lower")
  }
}

loc <- as.data.frame(matrix(nrow = nrow(pca$x), ncol = 2))
for(i in 1:nrow(loc)){
  loc[i,1] <- find.loc(pca$x[i,1:2][[1]],pca$x[i,1:2][[2]])
  loc[i,2] <- model.assignment[[i]]
}




ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = loc$V1)) +
  geom_point(size=2) + 
  labs(x="PC1 (36.3%)", y="PC2 (12.3%)", col="Treatment") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    legend.position="none") +
  geom_segment(aes(x = line.1.df[1,1], y = line.1.df[1,2], xend = line.1.df[nrow(line.1.df),1], yend = line.1.df[nrow(line.1.df),2] )) +
  geom_segment(aes(x = line.2.df[1,1], y = line.2.df[1,2], xend = line.2.df[nrow(line.2.df),1], yend = line.2.df[nrow(line.2.df),2] ))+
  scale_colour_manual(values=cbbPalette)


loc.index <- which(loc[,1] =="mid" )



ggplot(as.data.frame(pca$x[loc.index,]), aes(x = PC1, y = PC2, color = model.assignment[loc.index])) +
  geom_point(size=2) + 
  labs(x="PC1 (36.3%)", y="PC2 (12.3%)", col="Treatment") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    legend.position="none")+
  scale_colour_manual(values=cbbPalette)

##### make a linear model of subsetted data

#lmpca <- lm(PC2~PC1, data =as.data.frame(pca$x[loc.index,]) )
#summary(lmpca)

ggplot(as.data.frame(pca$x[loc.index,]), aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_abline(intercept = -5.5, slope = 2, col = "red") +
  theme_classic()+
  scale_colour_manual(values=cbbPalette)


####### project points to regression line
library(LearnGeom)


mk.line(m = 2, x = (-3):9, b = -5.5)

summary(lmpca)
line.start <- c(-3,-11.5)
line.end <-c(9,12.5)

line <- CreateLinePoints(line.start,line.end)



projected <- as.data.frame(matrix(nrow = nrow(pca$x[loc.index,1:2]), ncol = 3))

colnames(projected) <- c("X", "Y", "Assignment")

for(i in 1:nrow(pca$x[loc.index,1:2])){
  p1 <- pca$x[loc.index,1:2][i,]
  projected[i,1] <- ProjectPoint(p1, line)[[1]]
  projected[i,2] <- ProjectPoint(p1, line)[[2]]
  projected[i,3] <- model.assignment[loc.index][i]
}

#View(projected)

ggplot(as.data.frame(pca$x[loc.index,1:2]), aes(x = PC1, y = PC2))+geom_point()

ggplot(projected, aes(x = X, y = Y, color = model.assignment[loc.index])) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+
  scale_colour_manual(values=cbbPalette)

projected.IFNg <- projected[projected$Assignment == "IFNg",]
projected.IL4 <- projected[projected$Assignment == "IL4",]
projected.IFNg_IL4 <- projected[projected$Assignment == "IFNg_IL4",]


#hist(projected.IFNg[,1], col = "red")
#hist(projected.IL4[,1], col = "blue")
#hist(projected.IFNg_IL4[,1], col = "green")
#hist(projected.IL4[,1], col = "blue", add = T)
#hist(projected.IFNg_IL4[,1], col = "green", add = T)

projected.all <- rbind(projected.IFNg,projected.IL4,projected.IFNg_IL4 )

ggplot(projected.all, aes(x=X, fill=Assignment)) +
  geom_histogram( color='#e9ecef', alpha=0.6, position='identity', bins = 15) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))  +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"))+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"))

ggplot(projected.all, aes(x=X, fill=Assignment)) +
  geom_density( color='#e9ecef', alpha=0.6, position='identity') +
  xlab("PC1") +
  ylab("Density") + 
  theme_classic() +
  theme( 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))

######################################################################################################v
write.table(circ, file = "/Users/c-clausb/Desktop/Netact/object/circ.txt", quote = F, sep = "\t", row.names = F)








