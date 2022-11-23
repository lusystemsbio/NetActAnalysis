
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


inferLinks <- function(xdata=mydata, regdb, miTh = 0.1, nbins = 8,  method = "spearman"){
  # regulators=names(regdb)
  # xdata=mydata
  # miTh = 0.1; nbins = 3; method = "spearman"
  # 
  require(infotheo)
  require(reshape2)
  
  # retain expression data only for relevant TFs and their targets, whose expressions are available 
  dim(xdata)
  xdata <- xdata[rownames(xdata) %in% union(names(regdb), unlist(regdb)), ]
  dim(xdata)
  
  # Calculate the MI between pairs of nodes (retained TFs and targets) from the expression data
  #miMat = discretize(t(actMat), disc="equalfreq", nbins = 8) 
  miMat = discretize(t(xdata), disc="equalfreq", nbins = nbins)  
  miMat = mutinformation(miMat, method = "shrink") # calculate MI between node pairs
  diag(miMat) = 0
  actLinks = melt(miMat) 
  corMat = cor(t(xdata), method = method) # calculate correlation between nodes
  corLinks = melt(corMat) 
  
  # Retain the activation links whose sources are TFs in the regulons
  actLinks <- actLinks[actLinks$Var1 %in% names(regdb), ]
  dim(actLinks)
  # Retain the activation links whose targets are targets in regulons only
  actLinks <- actLinks[actLinks$Var1 %in% unlist(regdb), ]
  dim(actLinks)
  
  
  # construct the links by subsetting TF-target relations 
  tf_links <- actLinks[actLinks$value > miTh , c(1,2,3)]
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
  # dim(tf_links) 
  # 
  # return(tf_links)
  
  
  tf_links.dpi <- applyDPI(tf_links, miMat) 
  dim(tf_links.dpi) 
  tf_links.list <- list()
  tf_links.list[['nodpi']]  <- tf_links
  tf_links.list[['dpi']] <- tf_links.dpi
    
  return(tf_links.list)
}



# This function calculates precision (positive predictive value) and recall 
# for a list of MI TSH 
calPrecisionRecall <- function(inferredRel.df, regNet.df, miTSH.list){
  ppvRcall <- as.data.frame(matrix(nrow = length(miTSH.list), ncol = 2))
  rownames(ppvRcall) <- miTSH.list
  colnames(ppvRcall) <- c('ppv', 'rCall')
  
  for(i in 1:length(miTSH.list)) {
    inferredNet.df = inferredRel.df[inferredRel.df$mi>=miTSH.list[i], c(1, 2, 4)]
    dim(inferredNet.df)
    dupStatus <- duplicated(inferredNet.df )
    sum(dupStatus)
    
    inferredNet.df <- format_source_target_nodes(inputNet.df=inferredNet.df)
    dupStatus <- duplicated(inferredNet.df )
    sum(dupStatus)
    dim(inferredNet.df)
    
    inferredNet.retained.df <- inferredNet.df[, c(1, 2)] 
    colnames(inferredNet.retained.df) <- colnames(regNet.df) 
    dupStatus <- duplicated(inferredNet.retained.df)
    inferredNet.retained.df <- inferredNet.retained.df[!dupStatus,] # remove duplicated interactions
    
    combInt.df <- rbind(regNet.df, inferredNet.retained.df) 
    dupStatus <- duplicated(combInt.df)
    
    if(sum(dupStatus)>0) {
        rCall <- sum(dupStatus)/nrow(regNet.df)   
        ppv <- sum(dupStatus)/nrow(inferredNet.df) 
    }else {
        rCall <- ppv <- NA 
    }
    ppvRcall[i,] <- c(ppv, rCall)
  }  
  return(ppvRcall)
}

