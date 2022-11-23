
format_source_target_nodes <- function(inputNet.df){
  tmp.df <- as.data.frame(matrix(nrow = nrow(inputNet.df), ncol = ncol(inputNet.df)))
  for(i in 1:nrow(inputNet.df)){
    #i <- 1
    #print(i)
    myR <- inputNet.df[i, ] 
    if(as.character(myR[1]) > as.character(myR[2]) ) {
      tmp.df[i, ] <- c(as.character(myR[2]), as.character(myR[1]), as.numeric(myR[3]))
    }else tmp.df[i, ] <- myR
  }
  colnames(tmp.df) <- colnames(inputNet.df)
  return(tmp.df)
}

