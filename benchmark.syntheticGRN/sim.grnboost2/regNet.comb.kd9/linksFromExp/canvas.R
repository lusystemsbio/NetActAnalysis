

inputNet.df <- inferredNet.df

nodes.all <- unique(union(inputNet.df$TF, inputNet.df$target))
nodes.all <- sort(nodes.all)
length(nodes.all)


# tmp2.df <- format_source_target_nodes(inputNet.df = inferredNet.df) 
tmp3.df <- format_source_target_nodes(inputNet.df = regNet.df) 
 
regNet.df 

nodes.all.regNet <- sort(unique(union(regNet.df$SOURCE, regNet.df$TARGET)))
length(nodes.all.regNet) 

