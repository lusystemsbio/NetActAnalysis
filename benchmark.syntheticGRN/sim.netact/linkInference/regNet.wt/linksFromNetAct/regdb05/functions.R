"
  Swap source and target nodes to make them ordered lexicographically
"
format_source_target_nodes <- function(inputNet.df){
  #inputNet.df = regNet.df
  tmp.df <- as.data.frame(matrix(nrow = nrow(inputNet.df), ncol = ncol(inputNet.df)))
  for(i in 1:nrow(inputNet.df)){
    #i <- 1 
    myR <- inputNet.df[i, ] 
    if(as.character(myR[1]) > as.character(myR[2])) {
      tmp.df[i, ] <- c(as.character(unlist(myR[2])), as.character(unlist(myR[1])))
    }else tmp.df[i, ] <- c(as.character(unlist(myR[1])), as.character(unlist(myR[2])))
  }
  colnames(tmp.df) <- colnames(inputNet.df)
  return(tmp.df)
}


#' @export
#' @title Apply data processing inequality
#' @description Remove the interactions from a triangle which have lowest 
#' interaction score.
#' @param tfLinks Dataframe containing the interactions as source (character),
#'  target (character), type (integer).
#' @param miMat numeric matrix Interaction scores based on mutual information or 
#' correlation. 
#' @param miDiff numeric (0-1) Default 0.0 (optional) Minimum difference 
#' between mutual information of a triangle for the edge to be removed.
#' @param minMiTh numeric (0-1) Default 0.5. Minimum value of MI for an interaction 
#' which will not be removed. 
#' @return data.frame containing the filtered interactions.
#'   
applyDPI <- function(tfLinks = tfLinks, miMat = miMat, miDiff = 0, minMiTh = 0.5){
  print(miDiff)
  print(minMiTh)
  edgeKeep <- vector(mode = "integer", length = length(tfLinks[,1]))
  edgeRemove <- vector(mode = "integer", length = length(tfLinks[,1]))
  for(i in seq_along(tfLinks[,1])){
    srcGene <- tfLinks[i,1]
    tgtGene <- tfLinks[i,2]
    allTgts <- tfLinks[which(tfLinks[,1] == srcGene),2]
    for(j in seq_along(allTgts)){
      threeGenes <- c(srcGene, tgtGene, allTgts[j])
      if(length(unique(threeGenes)) == 3){
        thirdGenes <- union(tfLinks[which(tfLinks[,1] == allTgts[j]),2],
                            tfLinks[which(tfLinks[,2] == allTgts[j]),1])
        if((srcGene %in% thirdGenes) & (tgtGene %in% thirdGenes)){
          tmp <- miMat[c(srcGene,tgtGene,allTgts[j]),c(srcGene, 
                                                       tgtGene,allTgts[j])]
          minValue <- which.min(c(tmp[2,1],tmp[3,1],tmp[3,2]))
          minMi <- min(c(tmp[2,1],tmp[3,1],tmp[3,2]))
          maxMi <- max(c(tmp[2,1],tmp[3,1],tmp[3,2]))
          edge <- list()
          edge[[1]] <- which((tfLinks[,1] == srcGene) & (tfLinks[,2] == tgtGene) |
                               (tfLinks[,2] == srcGene) & (tfLinks[,1] == tgtGene)
          )
          edge[[2]] <- which((tfLinks[,1] == srcGene) & 
                               (tfLinks[,2] == allTgts[j]) | 
                               (tfLinks[,2] == srcGene) & 
                               (tfLinks[,1] == allTgts[j]))
          edge[[3]] <- which((tfLinks[,1] == tgtGene) & 
                               (tfLinks[,2] == allTgts[j]) | 
                               (tfLinks[,2] == tgtGene) & 
                               (tfLinks[,1] == allTgts[j]))
          
          #   edgeKeep[edge[1]] <- edgeKeep[edge[1]] + 1
          #    edgeKeep[edge[2]] <- edgeKeep[edge[2]] + 1
          #    edgeKeep[edge[3]] <- edgeKeep[edge[3]] + 1
          
          #   edgeKeep[edge[minValue]] <- edgeKeep[edge[minValue]] - 1
          #   edgeRemove[edge[minValue]] <- edgeRemove[edge[minValue]] - 1
          if((maxMi - minMi) > miDiff){
            if(minMi < minMiTh){
              edgeRemove[edge[[minValue]]] <- -1
            }
          }
        }
      }
    }
  }
  
  tfLinks <- tfLinks[-(which(edgeRemove < 0)),]
  return(tfLinks)
}


#' @export
#' @title Apply data processing inequality
#' @description Remove the interactions from a triangle which have lowest 
#' interaction score.
#' @param tfLinks Dataframe containing the interactions as source (character),
#'  target (character), type (integer).
#' @param miMat numeric matrix Interaction scores based on mutual information or 
#' correlation. 
#' @param miDiff numeric (0-1) Default 0.0 (optional) Minimum difference 
#' between mutual information of a triangle for the edge to be removed.
#' @param minMiTh numeric (0-1) Default 0.5. Minimum value of MI for an interaction 
#' which will not be removed. 
#' @return data.frame containing the filtered interactions.
#'   

applyDPI <- function(tfLinks = tfLinks, miMat = miMat, miDiff = 0, minMiTh = 0.5){
  print(miDiff)
  print(minMiTh)
  edgeKeep <- vector(mode = "integer", length = length(tfLinks[,1]))
  edgeRemove <- vector(mode = "integer", length = length(tfLinks[,1]))
  for(i in seq_along(tfLinks[,1])){
    srcGene <- tfLinks[i,1]
    tgtGene <- tfLinks[i,2]
    allTgts <- tfLinks[which(tfLinks[,1] == srcGene),2]
    for(j in seq_along(allTgts)){
      threeGenes <- c(srcGene, tgtGene, allTgts[j])
      if(length(unique(threeGenes)) == 3){
        thirdGenes <- union(tfLinks[which(tfLinks[,1] == allTgts[j]),2],
                            tfLinks[which(tfLinks[,2] == allTgts[j]),1])
        if((srcGene %in% thirdGenes) & (tgtGene %in% thirdGenes)){
          tmp <- miMat[c(srcGene,tgtGene,allTgts[j]),c(srcGene, 
                                                       tgtGene,allTgts[j])]
          minValue <- which.min(c(tmp[2,1],tmp[3,1],tmp[3,2]))
          minMi <- min(c(tmp[2,1],tmp[3,1],tmp[3,2]))
          maxMi <- max(c(tmp[2,1],tmp[3,1],tmp[3,2]))
          edge <- list()
          edge[[1]] <- which((tfLinks[,1] == srcGene) & (tfLinks[,2] == tgtGene) |
                               (tfLinks[,2] == srcGene) & (tfLinks[,1] == tgtGene)
          )
          edge[[2]] <- which((tfLinks[,1] == srcGene) & 
                               (tfLinks[,2] == allTgts[j]) | 
                               (tfLinks[,2] == srcGene) & 
                               (tfLinks[,1] == allTgts[j]))
          edge[[3]] <- which((tfLinks[,1] == tgtGene) & 
                               (tfLinks[,2] == allTgts[j]) | 
                               (tfLinks[,2] == tgtGene) & 
                               (tfLinks[,1] == allTgts[j]))
          
          #   edgeKeep[edge[1]] <- edgeKeep[edge[1]] + 1
          #    edgeKeep[edge[2]] <- edgeKeep[edge[2]] + 1
          #    edgeKeep[edge[3]] <- edgeKeep[edge[3]] + 1
          
          #   edgeKeep[edge[minValue]] <- edgeKeep[edge[minValue]] - 1
          #   edgeRemove[edge[minValue]] <- edgeRemove[edge[minValue]] - 1
          if((maxMi - minMi) > miDiff){
            if(minMi < minMiTh){
              edgeRemove[edge[[minValue]]] <- -1
            }
          }
        }
      }
    }
  }
  
  tfLinks <- tfLinks[-(which(edgeRemove < 0)),]
  return(tfLinks)
}


inferLinks <- function(xdata=mydata, regdb, miTh = 0.1, nbins = 3, 
                       method = "spearman", DPI = FALSE){
  require(infotheo)
  require(reshape2)
  
  # retain expression data only for relevant TFs and their targets, whose expressions are available 
  dim(xdata)
  sum(rownames(xdata) %in% union(names(regdb), unlist(regdb)))
  xdata <- xdata[rownames(xdata) %in% union(names(regdb), unlist(regdb)), ]
  dim(xdata)
  
  # Calculate the MI between pairs of nodes (retained TFs and targets) from the expression data
  #miMat = discretize(t(actMat), disc="equalfreq", nbins = 8) 
  miMat = discretize(t(xdata), disc="equalfreq", nbins = nbins)  
  miMat = mutinformation(miMat, method = "shrink") # calculate MI between node pairs
  diag(miMat) = 0 # set self interaction to 0
  actLinks = melt(miMat) 
  corMat = cor(t(xdata), method = method) # calculate correlation between nodes
  corLinks = melt(corMat) 
  
  # Retain the activation links whose sources are TFs in the regulons
  actLinks <- actLinks[actLinks$Var1 %in% names(regdb), ]

  dim(actLinks) 
  
  # Retain the activation links whose targets are target genes in regulons only
  # actLinks <- actLinks[actLinks$Var2 %in% unlist(regdb), ]
  # dim(actLinks)
  
  
  # construct the links by subsetting TF-target relations 
  #tf_links <- actLinks[actLinks$value > miTh , c(1,2,3)] 
  tf_links <- actLinks[actLinks$value >= miTh , c(1,2,3)] 
  colnames(tf_links) <- c('from', 'to', 'mi')
  dim(tf_links)
  
  
  # determine the directions based on correlations
  for (i in 1:nrow(tf_links)){ 
    source_name = tf_links$from[i]
    target_name = tf_links$to[i]
    tf_links$relation[i] <- ifelse(corMat[source_name, target_name] > 0, 1, 2)
  } 
  rownames(tf_links) = NULL 
  
  # tf_links <- applyDPI(tf_links, miMat)
  dim(tf_links) 
  if(DPI){
    tf_links <- applyDPI(tf_links, miMat)
  }
  dim(tf_links) 

  return(tf_links)
}


# This function calculates precision (positive predictive value) and recall 
# for a list of MI TSH 
calPrecisionRecall <- function(inferredRel.df, regNet.df, miTSH.list){
  ppvRcall <- as.data.frame(matrix(nrow = length(miTSH.list), ncol = 2))
  rownames(ppvRcall) <- miTSH.list
  colnames(ppvRcall) <- c('ppv', 'rCall')
  #i <- 9
  for(i in 1:length(miTSH.list)) {
    #print(i)
    inferredNet.df = inferredRel.df[inferredRel.df$mi>=miTSH.list[i], c(1, 2, 4)]
    if (nrow(inferredNet.df)==0){
      ppvRcall[i,] <- c(NA, NA)
      break
    }
    
    inferredNet.df <- inferredNet.df[, c(1, 2)]  # keep only source and target
    inferredNet.df <- format_source_target_nodes(inputNet.df=inferredNet.df)
    
    colnames(inferredNet.df) <- colnames(regNet.df) 
    dupStatus <- duplicated(inferredNet.df)
    inferredNet.retained.df <- inferredNet.df[!dupStatus,] # remove duplicated interactions
    
    combInt.df <- rbind(regNet.df, inferredNet.retained.df) 
    dupStatus <- duplicated(combInt.df)
    if(sum(dupStatus)>0){
      rCall <- sum(dupStatus)/nrow(regNet.df)   
      ppv <- sum(dupStatus)/nrow(inferredNet.retained.df) 
    }else rCall <- ppv <- NA
    ppvRcall[i,] <- c(ppv, rCall)
  }  
  return(ppvRcall)
}

# calculate confusion matrix - incomplete 
calConfusionMatrix <- function(inferredRel.df, regNet.df, miTSH.list){
  #inferredRel.df=inferredRel.df
  #regNet.df=regNet.df
  #miTSH.list=miTSH.list
  
  confusionMat <- as.data.frame(matrix(nrow = length(miTSH.list), ncol = 2))
  
  ppvRcall <- as.data.frame(matrix(nrow = length(miTSH.list), ncol = 2))
  rownames(ppvRcall) <- miTSH.list
  colnames(ppvRcall) <- c('ppv', 'rCall')
  i <- 3
  for(i in 1:length(miTSH.list)) {
    #print(i)
    inferredNet.df = inferredRel.df[inferredRel.df$mi>=miTSH.list[i], c(1, 2, 4)]
    #tmp.df <- inferredRel.df[inferredRel.df$mi<miTSH.list[i], c(1, 2, 4)]
    if (nrow(inferredNet.df)==0){
      ppvRcall[i,] <- c(NA, NA)
      break
    }
    
    inferredNet.df <- inferredNet.df[, c(1, 2)]  # keep only source and target
    inferredNet.df <- format_source_target_nodes(inputNet.df=inferredNet.df)
    
    colnames(inferredNet.df) <- colnames(regNet.df) 
    dupStatus <- duplicated(inferredNet.df)
    inferredNet.retained.df <- inferredNet.df[!dupStatus,] # remove duplicated interactions
    
    combInt.df <- rbind(regNet.df, inferredNet.retained.df) 
    dupStatus <- duplicated(combInt.df)
    if(sum(dupStatus)>0){
      rCall <- sum(dupStatus)/nrow(regNet.df)   
      ppv <- sum(dupStatus)/nrow(inferredNet.retained.df) 
    }else rCall <- ppv <- NA
    ppvRcall[i,] <- c(ppv, rCall)
  }  
  return(ppvRcall)
}




# Specific functions to calculate stats for new links 
#===========================================
regdb2String <- function(regdb.cur){
  regdbList.cur <- sapply(names(regdb.cur), function(tf){
    s <- sapply(regdb.cur[[tf]], function(tg) paste(tf, tg ,sep = ".", collapse = NULL))
    return(as.character(s))
  } 
  ) 
  return(regdbList.cur)
}

pairs2regdbFormat <- function(tfTargetPairs){
  tmp <- list()
  for(s in tfTargetPairs){
    s1 <- strsplit(s, split = '.', fixed = T)[[1]] 
    if(is.null(tmp[[s1[1]]])) tmp[[s1[1]]] <- s1[2] 
    else tmp[[s1[1]]] <- c(tmp[[s1[1]]], s1[2])  
  }
  return(tmp)
}

# construction regulatory network from regulon DB 
construct_regnet_from_regdb <- function(regdb){
  regNet <- NULL 
  for(tf in names(regdb)){
    targets.cur <- regdb[[tf]]
    tmp.df <- as.data.frame(matrix(nrow = length(targets.cur), ncol = 2)) 
    tmp.df[,1] <- rep(tf, length(targets.cur))
    tmp.df[,2] <- targets.cur
    regNet <- rbind(regNet, tmp.df)
  }
  colnames(regNet) <- c('SOURCE', 'TARGET') 
  return(regNet)
}

# Keep only thgose links in the inferred links that are NOt in the perturbed regulons
create_new_predictions <- function(inferredRel.df, regNet.pert){
  inferredRel.tmp <- inferredRel.df[,c(1,2)]
  dim(inferredRel.tmp)
  
  colnames(inferredRel.tmp) <- colnames(regNet.pert)
  combInt.df <- rbind(regNet.pert, inferredRel.tmp) 
  
  # Find the links in inferred link set that are in the perturbed regulon
  dupStatus <- duplicated(combInt.df) 
  
  # Find the links in the inferred link set that are not in perturbed regNet:
  dupStatus.n <- !dupStatus[(nrow(regNet.pert)+1):nrow(combInt.df)]
  new_predictions.df <- inferredRel.df[dupStatus.n,]
  return(new_predictions.df) 
}




df2Set <- function(myDF) {
  return(sapply(1:nrow(myDF), function(i) paste(myDF[i,1], myDF[i,2] ,sep = ".", collapse = NULL)))
}



# This function calculates precision (positive predictive value) and recall 
# for a list of MI TSH 
evaluatePredictions <- function(inferredRel.df, gtSet, miTSH.list){
  #inferredRel.df=inferredRel.df
  #gtSet <- gtSet
  #miTSH.list=miTSH.list
  
  gtSet.pos <- gtSet$regdbSet.pos
  gtSet.neg <- gtSet$regdbSet.neg
  
  #gtSet <- df2Set(gtLinks.df)   # convert the ground truth (GT) links to GT set
  
  stats.df <- as.data.frame(matrix(nrow = length(miTSH.list), ncol = 4))
  rownames(stats.df) <- miTSH.list
  colnames(stats.df) <- c('ppv', 'rCall', 'tpr',  'fpr') 
  #colnames(stats.df) <- c('ppv', 'rCall', 'Specificity', 'Sensitivity') 
  
  i <- 1
  for(i in 1:length(miTSH.list)) {
    #print(i)
    
    # Predicted links
    #----------------
    # find predicted links at the current TSH level: 
    inferredNet.df = inferredRel.df[inferredRel.df$mi>=miTSH.list[i], c(1, 2)]  
    if (nrow(inferredNet.df)==0){
      stats.df[i,] <- c(NA, NA, NA, NA)
      next
    }
    
    # keep only source and target for predicted links
    inferredNet.df <- format_source_target_nodes(inputNet.df=inferredNet.df) 
    dupStatus <- duplicated(inferredNet.df)
    inferredNet.retained.df <- inferredNet.df[!dupStatus,] # remove duplicated interactions
    inferredLinkSet <- df2Set(inferredNet.retained.df) # convert the predicted links to set of pairs (source, target)
    
    # True negatives 
    #---------------
    # calculate precision and recall using the inferred links 
    # precision/positive predictive value = TP/(TP+FP) ==> TP/(predictions) 
    # recall = TP/(TP+FN) ==> TP/P 
    
    # True positive rate or sensitivity is the same as recall
    # False positive rate or specificity: fpr = FP/N ==> FP/(FP+TN)
    tpSet <- intersect(gtSet.pos, inferredLinkSet) # true positive set
    fpSet <- setdiff(inferredLinkSet, gtSet.pos)# false positive set 
    # tnSet or true negative set is calculated before
    
    
    if(length(tpSet)>0){ # if true positive is not empty 
      rCall <- length(tpSet)/length(gtSet.pos) # gtSet: truePositiveSet + falseNegativeSet  
      ppv <- length(tpSet)/length(inferredLinkSet) # inferredLinkSet: predicted set
      
      # calculate roc values: TPR, FPR 
      # TPR = TP/P ==> TP/(TP+FN)
      # FPR = TP/N ==> TP/(FP+TN) 
      
      # calculate TPR
      tpRate <- rCall # TPR/Sensitivity is the same as recall 
      
      # calculate FPR i.e. 1-Specificity
      fpRate <- length(fpSet)/length(gtSet.neg) # FPR/
      #tnRate <- 1-fpRate # Specificity 
    }else{
      rCall <- ppv <- NA  
      tpRate  <- fpRate <- NA
    } 
    #stats.df[i,] <- c(ppv, rCall, tnRate, tpRate) 
    stats.df[i,] <- c(ppv, rCall, tpRate, fpRate) 
  } 
  return(stats.df)
}



