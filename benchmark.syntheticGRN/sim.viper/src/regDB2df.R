
# Convert toymodel regulon DB to doRothEA format 
regDB2df <- function(regdb){
  # convert the regulon db to df (doRothEA format)
  targetCount <- sapply(names(regdb), function(tf) length(regdb[[tf]]))
  noTargets <- sum(targetCount)
  
  regulons <- as.data.frame(matrix(nrow = length(names(regdb)), ncol =  4))
  colnames(regulons) <- c('tf', 'confidence', 'target', 'mor')
  
  cnt <- 1
  for(tf in names(regdb)){
    for (tg in regdb[[tf]]) {
      regulons[cnt, ] <- c(tf, 'C', tg, 1) 
      cnt <- cnt + 1
    }
  }
  class(regulons$mor)
  regulons$mor <- as.numeric(regulons$mor)  
  return(regulons)
}


#==================================================#
obtain_targetStat <- function(targets, vectorSize){
  targetStatus <- vector(mode = 'character', length = vectorSize)
  count <- 1
  for(tg in targets){ 
    targetStatus[count:(count+1)] <- c(tg, 1) 
    #targetStatus[count:(count+1)] <- c(tg, 100) 
    count <- count + 2
  }
  while (count<vectorSize) {
    #targetStatus[count:(count+1)] <- c('', NA)
    targetStatus[count:(count+1)] <- c(" ", 0.9 )
    count <- count + 2
  } 
  return(targetStatus)
}



#' MetaVIPER implementation that will perform a weighted stouffer integration based on highest NES.
#' 
#' @param ges Gene Expression Signature (features X samples)
#' @param net.list List object with the networks to be used
#' @param use.nets Optional argument to sslect the top n networks. If not specified, all networks are used.  
#' @param ret.weights Optional argument to return the network weight matrix as well as the VIPER matrix. FALSE by default.
#' @return Either a viper matrix, or a list with a viper matrix and the network weight matrix.
WeightedVIPER <- function(ges, net.list, use.nets, ret.weights = FALSE) {
  require(viper)
  num.nets <- length(net.list)
  num.samps <- ncol(ges)
  ## create weight matrix
  w.mat <- matrix(0L, nrow = num.nets, ncol = ncol(ges))
  colnames(w.mat) <- colnames(ges); rownames(w.mat) <- names(net.list)
  ## run VIPER with each network
  print('Generating VIPER matrices...')
  vip.list <- list()
  for (i in 1:num.nets) {
    vip.list[[i]] <- viper(ges, net.list[i], method = 'none')
  }
  names(vip.list) <- names(net.list)
  ## count for each gene
  print('Generating weights...')
  uni.genes <- unique(unlist(lapply(vip.list, rownames)))
  for (g in uni.genes) {
    for (s in 1:num.samps) {
      nes.vals <- unlist(lapply(vip.list, function(x){ 
        if (g %in% rownames(x)) {
          return(x[g,s])
        } else {
          return(0)
        }}))
      max.ind <- which.max(abs(nes.vals))
      w.mat[max.ind, s] <- w.mat[max.ind, s] + 1
    }
  }
  ## integration
  print('Integrating...')
  int.mat <- matrix(0L, nrow = length(uni.genes), ncol = num.samps)
  rownames(int.mat) <- uni.genes; colnames(int.mat) <- colnames(ges)
  for (g in uni.genes) {
    for (s in 1:num.samps) {
      nes.vals <- unlist(lapply(vip.list, function(x){ 
        if (g %in% rownames(x)) {
          return(x[g,s])
        } else {
          return(NA)
        }}))
      w.vals <- w.mat[,s][!is.na(nes.vals)]
      w.vals <- w.vals / sum(w.vals)
      nes.vals <- nes.vals[!is.na(nes.vals)]
      # if use.nets are specified, subset to the top n (or use all if n > length)
      if (!missing(use.nets)) {
        w.order <- order(w.vals, decreasing = TRUE)
        w.vals <- w.vals[ w.order[1:min(length(w.order), use.nets)] ]
        nes.vals <- nes.vals[ w.order[1:min(length(nes.vals), use.nets)] ]
      }
      int.mat[g,s] <- sum(nes.vals * w.vals) / sqrt(sum(w.vals**2))
    }
  }
  ## return
  if (ret.weights) {
    return( list('viper' = int.mat, 'weights' = w.mat) )
  } else {
    return( int.mat )
  }
} 


