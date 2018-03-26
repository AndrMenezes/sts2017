# Função de sobrevivência/distribuição conjunta ---------------------------

pHBE <- function(x, y, mu1, mu2, r, lower.tail = FALSE)
{
  Sxy <- exp(-((x / mu1) ^ r + (y / mu2) ^ r) ^ (1 / r))
  if(lower.tail) return(Sxy) else return(1 - Sxy)
}


# Função densidade conjunta -----------------------------------------------

dHBE <- function(x, y, mu1, mu2, r)
{
  fxy <- (x * y) ^ (r - 1) / ((mu1 * mu2) ^ r) * ((x / mu1) ^ r + (y / mu2) ^ r) ^ (1 / r - 2) *
    (r - 1 + ((x / mu1) ^ r + (y / mu2) ^ r) ^ (1 / r)) * pHBE(x, y, mu1, mu2, r, T)
  return(fxy)
}


# Função densidade condicional de Y dado X --------------------------------

dcondyHBE <- function(x, y, mu1, mu2, r)
{
  fy_x <- dHBE(x, y, mu1, mu2, r) / dexp(x, rate = 1 / mu1)
  return(fy_x)
}

pcondyHBE <- function(q, x, mu1, mu2, r)
{
  -0.1e1 / mu2 / x * ((x ^ r * mu2 ^ r + q ^ r * mu1 ^ r) ^ (-(-1 + 2 * r) / r) 
  * exp(-(-mu2 * x + (x ^ r * mu2 ^ r + q ^ r * mu1 ^ r) ^ (1 / r)) / mu1 / mu2)
  * q ^ r * mu2 ^ r * mu1 ^ r * x ^ r + (x ^ r * mu2 ^ r + q ^ r * mu1 ^ r) ^ (-(-1 + 2 * r) / r)
  * exp(-(-mu2 * x + (x ^ r * mu2 ^ r + q ^ r * mu1 ^ r) ^ (1 / r)) / mu1 / mu2) * mu2 ^ (2 * r)
  * x ^ (2 * r) - 0.1e1 / (x ^ r * mu2 ^ r + q ^ r * mu1 ^ r) * q ^ r * mu1 ^ r * mu2 * x
  - 0.1e1 / (x ^ r * mu2 ^ r + q ^ r * mu1 ^ r) * mu2 ^ (1 + r) * x ^ (1 + r))
}

# Quantil condicional de Y dado X -----------------------------------------

qcondyHBE <- function(p, x, mu1, mu2, r)
{
  Qy_x <- c()
  for(i in 1:length(x))
  {
    xi        <- x[i]
    pi        <- p[i]
    Fy_x      <- function(q) pcondyHBE(q, xi, mu1, mu2, r) - pi
    Qy_x[i]   <- uniroot(f = Fy_x, interval = c(0, 1000))$root
  }
  return(Qy_x)
}

# Geração usando método da condicional ------------------------------------

rHBE <- function(n, mu1, mu2, r)
{
  x <- rexp(n, 1 / mu1)
  u <- runif(n)
  y <- qcondyHBE(p = u, x = x, mu1 = mu1, mu2 = mu2, r = r)
  return(cbind(x = x, y = y))
}

# Função Log-Verossimilhança ----------------------------------------------
llike_HBE <- function(parms, x, y)
{
  mu1 <- parms[1]
  mu2 <- parms[2]
  r   <- parms[3] 
  
  ll <- sum(log(dHBE(x, y, mu1, mu2, r)))
  return(-ll)
}


# Coeficiente de correlação -----------------------------------------------
corr_HBE <- function(r)
{
  rho <- gamma(1 / r)^2 / (r * gamma(2 / r)) - 1
  return(rho)
}
# r = 10
# corr_HBE(r)
# dados <- rHBE(1000, 3, 2, r)
# cor(dados[, 1], dados[, 2])


# Confiabilidade ----------------------------------------------------------
R_HBE <- function(mles)
{
  mu1 <- mles[1]
  mu2 <- mles[2]
  r   <- mles[3] 
  R   <- mu2 ^ r / (mu1 ^ r + mu2 ^ r);
  return(R)
}




# Exemplo 1 ---------------------------------------------------------------

# library(maxLik)
# set.seed(666)
# dados <- rHBE(1000, 2, 3, 2)
# x     <- dados[, 1]
# y     <- dados[, 2]
# out   <- optim(par = c(2, 3, 2), fn = llike_HBE, method = "BFGS", x = x, y = y)
# out

# Exemplo 2 ---------------------------------------------------------------
# set.seed(666)
# dados <- rHBE(n = 100, mu1 = 5, mu2 = 10, r = 8)
# x     <- dados[, 1]
# y     <- dados[, 2]
# out   <- optim(par = c(5, 10, 8), fn = llike_HBE, method = "BFGS", x = x, y = y)
# out
# 
# write.table(x = dados, file = "exemplo-HBE.csv", quote = F, sep = ',', row.names = F)



