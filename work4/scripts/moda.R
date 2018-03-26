moda <- function(dados)
{
  ht=hist(dados, plot = F)
  mcl <- which.max(ht$counts) # i 
  li <- ht$breaks[mcl] # li
  width <- diff(ht$breaks[mcl+0:1]) # h
  counts <- c(0,ht$counts,0)
  delta <- abs(diff(counts[1+mcl+(-1:1)])) # denominador
  moda <- li+width*delta[1]/sum(delta)
  cols <- rep(5, length(ht$counts))
  cols[mcl] <- 3
  return(moda=moda)
}
