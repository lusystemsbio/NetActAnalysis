
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
      tmp.df[i, ] <- c(as.character(unlist(myR[2])), as.character(unlist(myR[1])), as.character(unlist(myR[3])))
    }else tmp.df[i, ] <- c(as.character(unlist(myR[1])), as.character(unlist(myR[2])), as.character(unlist(myR[3])))
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

# inferLinks <- function(xdata=mydata, regdb, miTh = 1.4, nbins = 8, method = "spearman"){
inferLinks <- function(xdata=mydata, regdb, miTh = 0.1, nbins = 8,  
                       method = "spearman", DPI = FALSE){
  #regulators=names(regdb)
  #xdata=mydata
  #miTh = 0.0; nbins = 3; method = "spearman"
  
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
  actLinks <- actLinks[actLinks$Var2 %in% unlist(regdb), ]
  dim(actLinks)
  
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
  #i <- 6
  for(i in 1:length(miTSH.list)) {
    inferredNet.df = inferredRel.df[inferredRel.df$mi>=miTSH.list[i], c(1, 2, 4)]
    if (nrow(inferredNet.df)==0){
      ppvRcall[i,] <- c(0, 0)
      break
    }
    inferredNet.df <- format_source_target_nodes(inputNet.df=inferredNet.df)
    
    inferredNet.retained.df <- inferredNet.df[, c(1, 2)] 
    colnames(inferredNet.retained.df) <- colnames(regNet.df) 
    dupStatus <- duplicated(inferredNet.retained.df)
    inferredNet.retained.df <- inferredNet.retained.df[!dupStatus,] # remove duplicated interactions
    
    combInt.df <- rbind(regNet.df, inferredNet.retained.df) 
    dupStatus <- duplicated(combInt.df)
    if(sum(dupStatus)>0) { 
        rCall <- sum(dupStatus)/nrow(regNet.df)   
        #ppv <- sum(dupStatus)/nrow(inferredNet.df) 
        ppv <- sum(dupStatus)/nrow(inferredNet.retained.df) 
        ppvRcall[i,] <- c(ppv, rCall)
    }else {
        ppvRcall[i,] <- c(NA, NA)
    }
  }  
  return(ppvRcall)
}






