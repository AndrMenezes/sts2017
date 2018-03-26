# Função de sobrevivência/distribuição conjunta ---------------------------

pbexp <- function(x, y, lambda1, lambda2, theta, lower.tail = FALSE)
{
  Sxy <- exp(- (lambda1 * x + lambda2 * y + theta * lambda1 * lambda2 * x * y))
  if(lower.tail) return(Sxy) else return(1 - Sxy)
}


# Função densidade conjunta -----------------------------------------------

dbexp <- function(x, y, lambda1, lambda2, theta)
{
  fxy <- ((1 - theta) * lambda1 * lambda2 + theta * lambda1^2 * lambda2 * x + theta * lambda1 * lambda2 ^ 2 * y + 
            theta ^ 2 * lambda1 ^ 2 * lambda2 ^ 2 * x * y) * pbexp(x, y, lambda1, lambda2, theta, TRUE)
  return(fxy)
}


# Função densidade condicional de Y dado X --------------------------------

dcondy <- function(x, y, lambda1, lambda2, theta)
{
  fy_x <- dbexp(x, y, lambda1, lambda2, theta) / dexp(x, rate = lambda1)
  return(fy_x)
}

# Quantil condicional de Y dado X -----------------------------------------

qcondy <- function(p, x, lambda1, lambda2, theta, exact = TRUE)
{
  Qy_x <- c()
  if(exact == FALSE)
  {
    for(i in 1:length(x))
    {
      xi        <- x[i]
      pi        <- p[i]
      integrand <- function(y) dbexp(xi, y, lambda1, lambda2, theta) / dexp(xi, rate = lambda1)
      Fy_x      <- function(q) integrate(integrand, lower = 0, upper = q)$value - pi
      Qy_x[i]   <- uniroot(f = Fy_x, interval = c(0, 100), maxiter = 10000)$root
    }
  }
  else
  {
    for(i in 1:length(x))
    {
      Qy_x[i]   <- -(lambda1 * theta * x[i] + theta * lamW::lambertWm1((-1 + p[i]) *
                    (lambda1 * theta * x[i] + 1) / theta * exp(-(lambda1 * theta * x[i] + 1) / theta)) + 1) / (lambda1 * theta * x[i] + 1) / lambda2 / theta
    }
  }
  return(Qy_x)
}

# Geração usando método da condicional ------------------------------------

rbexp <- function(n, lambda1, lambda2, theta)
{
  x <- rexp(n, lambda1)
  u <- runif(n)
  y <- qcondy(p = u, x = x, lambda1 = lambda1, lambda2 = lambda2, theta = theta, exact = TRUE)
  return(cbind(x = x, y = y))
}

# Função Log-Verossimilhança ----------------------------------------------
llike_bexp <- function(parms, x, y)
{
  lambda1 <- parms[1]
  lambda2 <- parms[2]
  theta   <- parms[3] 
  
  ll <- sum(log(dbexp(x, y, lambda1, lambda2, theta)))
  return(-ll)
}


# Coeficiente de correlação -----------------------------------------------
E1 <- function(x)
{
  A <- log((0.56146 / x + 0.65) * (1 + x))
  B <- x ^ 4 * exp(7.7 * x) * (2 + x) ^ (3.7)
  e <- (A ^ (-7.7) + B) ^ (-0.13)   
  return(e)
}
corr_bexp <- function(theta, exact = TRUE)
{
  if(exact) Ei  <- expint::expint_E1(1 / theta)
  else Ei <- E1(1 / theta)
  rho <- 1 / theta * exp(1 / theta) * Ei - 1
  return(rho)
}
theta = 0.8254
corr_bexp(theta)
corr_bexp(theta, exact = F)
dados <- rbexp(10000, 2, 3, theta)
cor(dados[, 1], dados[, 2])


# Confiabilidade ----------------------------------------------------------
R_bexp <- function(mles)
{
  lambda1 <- mles[1]
  lambda2 <- mles[2]
  theta   <- mles[3] 
  aux     <-	pnorm((lambda1 + lambda2) / (sqrt(2 * theta * lambda1 * lambda2)))
  R       <- 1 / 2 + sqrt(pi) * (lambda1 - lambda2) / (2 * sqrt(theta * lambda1 * lambda2)) *
             exp((lambda1 + lambda2)**2 / (4 * theta * lambda1 * lambda2)) * (1 - aux)
  return(R)
}




# Exemplo 1 ---------------------------------------------------------------

# library(maxLik)
# set.seed(666)
# dados <- rbexp(100, 2, 3, 0.5)
# x     <- dados[, 1]
# y     <- dados[, 2]
# out   <- optim(par = c(2, 3, 0.5), fn = llike_bexp, method = "BFGS")
# out

# Exemplo 2 ---------------------------------------------------------------
# set.seed(666)
# dados <- rbexp(n = 100, lambda1 = 2, lambda2 = 3, theta = 0.8)
# x     <- dados[, 1]
# y     <- dados[, 2]
# out   <- optim(par = c(2, 3, 0.8), fn = llike_bexp, method = "BFGS", x = x, y = y)
# out

# write.table(x = dados, file = "exemplo2.csv", quote = F, sep = ',', row.names = F)




# E1 <- function(x)
# {
#   A <- log((0.56146 / x + 0.65) * (1 + x))
#   B <- x ^ 4 * exp(7.7 * x) * (2 + x) ^ (3.7)
#   e <- (A ^ (-7.7) + B) ^ (-0.13)   
#   return(e)
# }
# expint_E1(theta)
# E1(theta)
# expint_En(theta, 1)




