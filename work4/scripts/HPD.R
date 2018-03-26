hpd.Chen <- function(x, alpha = 0.05)
{
  x    <- sort(x)
  n    <- length(x)
  namp <- (1 - alpha) * n
  amp  <- sapply(1:(n - namp), function(j)   diff(c(x[j], x[j + namp])))
  id   <- order(amp)[1]
  c(lower = x[id], upper = x[id + namp])
}