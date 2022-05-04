########### apply sign switching to macrophage data ##############


rm(list = ls())
setwd("/Users/c-clausb/Desktop/NetAct/reneedsomehelpfornetactmodeling")
# Using saved datasets - macrophageData.RDS, counts.rds, info.txt, 2019-03-07_circuit_file11c532b7432da.RDS
library("AnnotationDbi")
library("org.Mm.eg.db")
library(NetAct)
# columns(org.Mm.eg.db)
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
DErslt = MultiRNAseqDegs_limma(count_data, phenoData, compare_list)


library(sva)
e = DErslt$Overall$e
rownames(phenoData) = colnames(e)

# design = model.matrix(~1, data = pData(phenoData))

design <- model.matrix(~grouping)

combat_edata = ComBat(dat = as.matrix(e), batch = batch, mod = design, par.prior=T)


neweset = ExpressionSet(assayData = as.matrix(combat_edata), 
                        phenoData = phenoData)


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
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

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



#gsearslt_1 = TF_GSEA(mgs, DErslt$`UT-IL4_2h`, minSize=5, nperm = 10000, qval = T)
#gsearslt_2 = TF_GSEA(mgs, DErslt$`UT-IFNg_2h`, minSize=5, nperm = 10000, qval = T)
#gsearslt_3 = TF_GSEA(mgs, DErslt$`UT-IFNg_IL4_2h`, minSize=5, nperm = 10000, qval = T)
#gsearslt_4 = TF_GSEA(mgs, DErslt$`UT-IL4_4h`, minSize=5, nperm = 10000, qval = T)
#gsearslt_5 = TF_GSEA(mgs, DErslt$`UT-IFNg_4h`, minSize=5, nperm = 10000, qval = T)
#gsearslt_6 = TF_GSEA(mgs, DErslt$`UT-IFNg_IL4_4h`, minSize=5, nperm = 10000, qval = T)

#gsearslts <- list(gsearslt_1,gsearslt_2,gsearslt_3, gsearslt_4,gsearslt_5,gsearslt_6)
#names(gsearslts) <- c("gsearslt_1","gsearslt_2","gsearslt_3", "gsearslt_4","gsearslt_5","gsearslt_6")

#macrophageData <- list(neweset, DErslt, combat_edata,count_data, design, data_macro, gsearslts)
#names(macrophageData) <- c("neweset", "DErslt", "combat_edata","count_data", "design", "data_macro", "gsearslts")

#saveRDS(macrophageData, file = "macrophageData.RDS")

########### netact #################
macrophageData <- readRDS(file = "macrophageData.RDS")

qValue.threshold <- 0.01
######## for Qval 0.05 only #####
#qValue.threshold <- 0.05
######### adjusting qvalue of just one comparison #############
#tfs_1 = as.character(macrophageData$gsearslts$gsearslt_1$tf[macrophageData$gsearslts$gsearslt_1$qvals 
#                                                            < qValue.threshold])

tfs_1 = as.character(macrophageData$gsearslts$gsearslt_1$tf[macrophageData$gsearslts$gsearslt_1$qvals 
                                                            < .05])
##################


tfs_2 = as.character(macrophageData$gsearslts$gsearslt_2$tf[macrophageData$gsearslts$gsearslt_2$qvals 
                                                            < qValue.threshold])
tfs_5 = as.character(macrophageData$gsearslts$gsearslt_2$tf[macrophageData$gsearslts$gsearslt_5$qvals 
                                                            < qValue.threshold])
tfs_3 = as.character(macrophageData$gsearslts$gsearslt_3$tf[macrophageData$gsearslts$gsearslt_3$qvals 
                                                            < qValue.threshold])
tfs_4 = as.character(macrophageData$gsearslts$gsearslt_1$tf[macrophageData$gsearslts$gsearslt_4$qvals 
                                                            < qValue.threshold])
tfs_6 = as.character(macrophageData$gsearslts$gsearslt_3$tf[macrophageData$gsearslts$gsearslt_6$qvals 
                                                            < qValue.threshold])
tfs = sort(unique(as.character(c(tfs_1,tfs_2,tfs_3, tfs_4,tfs_5, tfs_6))))
#lengths(list(tfs_1,tfs_2,tfs_3, tfs_4,tfs_5, tfs_6,tfs))
#tfs

compare_list = c("UT-IL4_2h", "UT-IFNg_2h", "UT-IFNg_IL4_2h","UT-IL4_4h", "UT-IFNg_4h", "UT-IFNg_IL4_4h")
# acts_mat = TF_Activity(tfs, GSDB, counts, DErslt, with_weight = TRUE, useDatabaseSign = F, useCorSign = T, if_module = T)
act.me <- TF_Activity(tfs, mgs, macrophageData$neweset, macrophageData$DErslt$Overall, with_weight = TRUE, useDatabaseSign = F, useCorSign = T, if_module = F)
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
tfLinks1 = TF_Filter(actsMat, mgs, miTh = miTh, nbins =3, maxTf = 150, 
                     maxInteractions = 500,  corMethod = "s", 
                     useCor = T, removeSignalling = FALSE, DPI = F, 
                     miDiff = 0.0, minMiTh = minMiTh)
dim(tfLinks1)
tfs_network <- sort(union(tfLinks1[,1],tfLinks1[,2]))
tfs_network
#plot_network_v(tfLinks1)

tfsSelected <- rownames(acts_mat)[which(rownames(acts_mat) %in% sort(unique(as.character(c(tfs_2,tfs_5)))))]
actsMat <- acts_mat[tfsSelected,c("UT_R1", "UT_R2", "UT_R3",  "IFNg_2h_R1",  "IFNg_2h_R2",  "IFNg_2h_R3", "IFNg_4h_R1", "IFNg_4h_R2", "IFNg_4h_R3"   )]
tfLinks2 = TF_Filter(actsMat, mgs, miTh = miTh, nbins =3, maxTf = 150, 
                     maxInteractions = 500,  corMethod = "s", 
                     useCor = T, removeSignalling = FALSE, DPI = F, 
                     miDiff = 0.0, minMiTh = minMiTh)
dim(tfLinks2)
#plot_network_v(tfLinks2)

tfsSelected <- rownames(acts_mat)[which(rownames(acts_mat) %in% sort(unique(as.character(c(tfs_3,tfs_6)))))]
actsMat <- acts_mat[tfsSelected,c("UT_R1", "UT_R2", "UT_R3",  "IFNg_IL4_2h_R1",  "IFNg_IL4_2h_R2",  "IFNg_IL4_2h_R3", "IFNg_IL4_4h_R1", "IFNg_IL4_4h_R2", "IFNg_IL4_4h_R3"   )]
tfLinks3 = TF_Filter(actsMat, mgs, miTh = miTh, nbins =3, maxTf = 150, 
                     maxInteractions = 500,  corMethod = "s", 
                     useCor = T, removeSignalling = FALSE, DPI = F, 
                     miDiff = 0.0, minMiTh = minMiTh)
dim(tfLinks3)
# plot_network_v(tfLinks3)



tfLinks <- rbind(tfLinks1,tfLinks2, tfLinks3)
# plot_network_v(tfLinks)

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

plot_network_v(tfLinksN)

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
tfLinks1 = TF_Filter(actsMat, mgs, miTh = miTh, nbins =3, maxTf = 150, 
                     maxInteractions = 500,  corMethod = "s", 
                     useCor = T, removeSignalling = FALSE, DPI = F, 
                     miDiff = 0.0, minMiTh = minMiTh)
dim(tfLinks1)
tfs_network <- sort(union(tfLinks1[,1],tfLinks1[,2]))
tfs_network
#plot_network_v(tfLinks1)

tfsSelected <- rownames(acts_mat)[which(rownames(acts_mat) %in% sort(unique(as.character(c(tfs_2,tfs_5)))))]
actsMat <- acts_mat[tfsSelected,c("UT_R1", "UT_R2", "UT_R3",  "IFNg_2h_R1",  "IFNg_2h_R2",  "IFNg_2h_R3", "IFNg_4h_R1", "IFNg_4h_R2", "IFNg_4h_R3"   )]
tfLinks2 = TF_Filter(actsMat, mgs, miTh = miTh, nbins =3, maxTf = 150, 
                     maxInteractions = 500,  corMethod = "s", 
                     useCor = T, removeSignalling = FALSE, DPI = F, 
                     miDiff = 0.0, minMiTh = minMiTh)
dim(tfLinks2)
#plot_network_v(tfLinks2)

tfsSelected <- rownames(acts_mat)[which(rownames(acts_mat) %in% sort(unique(as.character(c(tfs_3,tfs_6)))))]
actsMat <- acts_mat[tfsSelected,c("UT_R1", "UT_R2", "UT_R3",  "IFNg_IL4_2h_R1",  "IFNg_IL4_2h_R2",  "IFNg_IL4_2h_R3", "IFNg_IL4_4h_R1", "IFNg_IL4_4h_R2", "IFNg_IL4_4h_R3"   )]
tfLinks3 = TF_Filter(actsMat, mgs, miTh = miTh, nbins =3, maxTf = 150, 
                     maxInteractions = 500,  corMethod = "s", 
                     useCor = T, removeSignalling = FALSE, DPI = F, 
                     miDiff = 0.0, minMiTh = minMiTh)
dim(tfLinks3)
# plot_network_v(tfLinks3)



tfLinks <- rbind(tfLinks1,tfLinks2, tfLinks3)
# plot_network_v(tfLinks)

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





plot_network_v(tfLinksN)
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

H1 =  ComplexHeatmap::Heatmap(row_norm(new_activity), col = circlize::colorRamp2(c(-2, 0, 2), c("green3", "white", "red")),
                              cluster_columns = F, cluster_rows = T, show_row_dend = FALSE, name = "Activity", column_title = "Activity",
                              row_names_gp = grid::gpar(fontsize = 10), column_names_gp = grid::gpar(fontsize = 7), clustering_method_rows = "ward.D2",
                              clustering_distance_rows = function(x) as.dist((1-cor(t(x), method = "spear"))/2)) 
# ,
#                              top_annotation = HeatmapAnnotation(data.frame(Treatment = substr(colnames(new_activity),1,(nchar(colnames(new_activity))-3))), col = list(Treatment = groupColors)))
gs = rownames(new_activity)
gc = colnames(new_activity)
H2 = ComplexHeatmap::Heatmap(row_norm(data[gs, gc]), col = circlize::colorRamp2(c(-2, 0, 2), c("green3", "white", "red")),
                             cluster_columns = F, cluster_rows = T, show_row_dend = FALSE,name = "Expression", column_title = "Expression",
                             row_names_gp = grid::gpar(fontsize = 10), column_names_gp = grid::gpar(fontsize = 7), clustering_method_rows = "ward.D2",
                             clustering_distance_rows = function(x) as.dist((1-cor(t(x), method = "spear"))/2))
#,
#                             top_annotation = HeatmapAnnotation(data.frame(Treatment = substr(colnames(new_activity),1,(nchar(colnames(new_activity))-3))), col = list(Treatment = groupColors)))
H = H1 + H2
draw(H, show_annotation_legend = FALSE)



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

ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = model.assignment)) +
  geom_point(size=5) + 
  labs(x="PC1 (36.6.6%)", y="PC2 (12.4%)", col="Treatment") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    legend.position="none"
  )








#ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = model.assignment)) + geom_point() + stat_ellipse()



#View(model.assignment)
model.assignment[model.assignment == "IFNg_IL4_2h" | model.assignment == "IFNg_IL4_4h"] <- "IFNg_IL4"
model.assignment[model.assignment == "IFNg_2h" | model.assignment == "IFNg_4h"] <- "IFNg"
model.assignment[model.assignment == "IL4_2h" | model.assignment == "IL4_4h"] <- "IL4"

#save(model.assignment, file = "/Users/c-clausb/Desktop/Netact/object/MacFinalModelAssignment")
load("/Users/c-clausb/Desktop/Netact/object/MacFinalModelAssignment")


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
  geom_segment(aes(x = line.2.df[1,1], y = line.2.df[1,2], xend = line.2.df[nrow(line.2.df),1], yend = line.2.df[nrow(line.2.df),2] ))

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
  geom_segment(aes(x = line.2.df[1,1], y = line.2.df[1,2], xend = line.2.df[nrow(line.2.df),1], yend = line.2.df[nrow(line.2.df),2] ))


loc.index <- which(loc[,1] =="mid" )



ggplot(as.data.frame(pca$x[loc.index,]), aes(x = PC1, y = PC2, color = model.assignment[loc.index])) +
  geom_point(size=2) + 
  labs(x="PC1 (36.3%)", y="PC2 (12.3%)", col="Treatment") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    legend.position="none")

##### make a linear model of subsetted data

#lmpca <- lm(PC2~PC1, data =as.data.frame(pca$x[loc.index,]) )
#summary(lmpca)

ggplot(as.data.frame(pca$x[loc.index,]), aes(x = PC1, y = PC2, color = model.assignment[loc.index])) +
  geom_point(size=2) + 
  labs(x="PC1 (36.3%)", y="PC2 (12.3%)", col="Treatment") +
  theme_bw() + theme(#panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    text = element_text(size=20),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    legend.position="none")+
  geom_abline(intercept = -5.5, slope = 2, col = "red")



ggplot(as.data.frame(pca$x[loc.index,]), aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_abline(intercept = -5.5, slope = 2, col = "red") +
  theme_classic()


####### project points to  line
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
        axis.title = element_text(size = 20))

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
        axis.title = element_text(size = 20))



ggplot(projected.all, aes(x=X, fill=Assignment)) +
  geom_density( color='#e9ecef', alpha=0.6, position='identity') +
  xlab("PC1") +
  ylab("Density") + 
  theme_classic() +
  theme( 
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20))

######################################################################################################v
#write.table(circ, file = "/Users/c-clausb/Desktop/Netact/object/circ.txt", quote = F, sep = "\t", row.names = F)









########## heatmap  #############

#rm(list = ls())

library(sRACIPE)

load("/Users/c-clausb/Desktop/Netact/object/gex_normal")
load("/Users/c-clausb/Desktop/Netact/object/MacFinalModelAssignment")


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










clustMethod = "ward.D2"
distType = "euclidean"
corMethod = "spearman"
assignedClusters = model.assignment
assayDataTmp <- gex_normal
col <-  grDevices::colorRampPalette(rev(
  RColorBrewer::brewer.pal(11, 'Spectral')))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
brc.colors <- gg_color_hue(4)
brc.colors <- brc.colors[c(2,3,1,4)]
col2 <- brc.colors





corCol <- stats::cor((assayDataTmp), method = corMethod)
distanceCol <- stats::as.dist((1 - corCol) / 2)
corRow <- stats::cor(t(assayDataTmp), method = corMethod)
distanceRow <- stats::as.dist((1 - corRow) / 2)
clustersCol <- stats::hclust(distanceCol, method = clustMethod)
ddCol <- as.dendrogram(clustersCol)
clustersRow <- stats::hclust(distanceRow, method = clustMethod)
ddRow <- stats::as.dendrogram(clustersRow)




clustNames <- unique(assignedClusters)
nClusters <- length(clustNames)
clustColors <- numeric(length(assignedClusters))
for(tmp1 in seq_len(length(clustColors))){
  clustColors[tmp1] <- which(clustNames == assignedClusters[tmp1] )
}
clustColors <- col2[clustColors]
names(clustColors) <- assignedClusters




tiff("/Users/c-clausb/Desktop/trash/mac_heat", width=2000,height=3000,res=400, compression="lzw")
gplots::heatmap.2((assayDataTmp),
                  Colv = ddCol,
                  Rowv = ddRow,
                  trace = "none",
                  col = col,
                  ColSideColors = clustColors,
                  cexRow = 0.6
)

dev.off()

