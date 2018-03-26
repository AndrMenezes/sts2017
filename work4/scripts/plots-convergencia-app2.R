# Análise de convergência das posteriori marginais ------------------------

rm(list = ls())
wd1 <- 'C:/Users/User/Dropbox/4° Série/Tópicos Especiais em Estatística/Trabalhos/Trabalho 4/Scripts' 
wd2 <- 'C:/Users/User/Dropbox/4° Série/Tópicos Especiais em Estatística/Trabalhos/Trabalho 4' 
source(paste0(wd1, '/' ,'MOEW.R'))
source(paste0(wd1, '/' ,'moda.R'))
FF  <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}

# Amostras simuladas versus a iteração ------------------------------------

it.plot <- function(par, npar, sex)
{
  n <- length(par)
  R <- range(par)
  #pdf(file = paste0(npar, '-st-', sex, '.pdf'), width = 8)
  jpeg(file = paste0(npar, '-st-', sex, '.jpeg'), quality = 100, res = 300, width = 8, height = 7, units = 'in')
  par(mar = c(3.2, 3.2, 1.0, 1.0), cex = 1.8)
  plot(par, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', type = 'l')
  axis(side = 1, at = seq(0, 20000, by = 5000), labels = FF(seq(0, 20000, by = 5000),0))
  axis(side = 2, at = seq(R[1], R[2], l = 5), labels = FF(seq(R[1], R[2], l = 5), 2))
  mtext(text = 'Iteração', side = 1, line = 2, cex = 1.8)
  if(npar == 'mu') mtext(text = expression(mu), side = 2, line = 2, cex = 1.8)
  if(npar == 'beta') mtext(text = expression(beta), side = 2, line = 2, cex = 1.8)
  if(npar == 'alpha') mtext(text = expression(alpha), side = 2, line = 2, cex = 1.8)
  graphics.off()
}

# Gráficos da densidade posteriori ----------------------------------------

densi.plot <- function(par, npar, sex, legend = F)
{
  densi <- density(par)
  Rx    <- range(densi$x)
  Ry    <- range(densi$y)
  pdf(file = paste0(npar, '-densi-', sex, '.pdf'), width = 8)
  # jpeg(file = paste0(npar, '-densi-', sex, '.jpeg'), width = 580, height = 480, quality = 100)
  par(mar = c(3.2, 3.2, 1.0, 1.0), cex = 1.8)
  plot(densi, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', type = 'l', main = '', lwd = 2)
  rug(par)
  axis(side = 1, at = seq(Rx[1], Rx[2], l = 5), labels = FF(seq(Rx[1], Rx[2], l = 5), 2))
  axis(side = 2, at = seq(Ry[1], Ry[2], l = 5), labels = FF(seq(Ry[1], Ry[2], l = 5), 2))
  mtext(text = 'Densidade', side = 2, line = 2, cex = 1.8)
  if(npar == 'mu') mtext(text = expression(mu), side = 1, line = 2, cex = 1.8)
  if(npar == 'beta')   mtext(text = expression(beta),  side = 1, line = 2, cex = 1.8)
  if(npar == 'alpha')  mtext(text = expression(alpha), side = 1, line = 2, cex = 1.8)
  abline(v = mean(par), col = 'red', lwd = 2)
  abline(v = median(par), col = 'blue', lwd = 2)
  abline(v = moda(par), col = 'green', lwd = 2)
  if(legend == TRUE)  legend('topright', legend = c('Média', 'Mediana', 'Moda'), col = c('red', 'blue', 'green'), bty = 'n', lwd = 2)
  graphics.off()
}


# MOEW - Female -----------------------------------------------------------
setwd(wd1)
posterioris <- read.delim('post2.txt')
head(posterioris)
setwd(wd2)

mus    <- posterioris[which(posterioris$gender == 'female'),]$mu
betas  <- posterioris[which(posterioris$gender == 'female'),]$beta
alphas <- posterioris[which(posterioris$gender == 'female'),]$alpha

it.plot(par = mus,    npar = 'mu',    sex = 'fem')
it.plot(par = betas,  npar = 'beta',  sex = 'fem')
it.plot(par = alphas, npar = 'alpha', sex = 'fem')

densi.plot(par = mus,    npar = 'mu',    legend = T, sex = 'fem')
densi.plot(par = betas,  npar = 'beta',  legend = T, sex = 'fem')
densi.plot(par = alphas, npar = 'alpha', legend = T, sex = 'fem')

# MOEW - Male -----------------------------------------------------------

mus    <- posterioris[which(posterioris$gender == 'male'),]$mu
betas  <- posterioris[which(posterioris$gender == 'male'),]$beta
alphas <- posterioris[which(posterioris$gender == 'male'),]$alpha

it.plot(par = mus,    npar = 'mu',    sex = 'man')
it.plot(par = betas,  npar = 'beta',  sex = 'man')
it.plot(par = alphas, npar = 'alpha', sex = 'man')

densi.plot(par = mus,    npar = 'mu',    legend = T, sex = 'man')
densi.plot(par = betas,  npar = 'beta',  legend = T, sex = 'man')
densi.plot(par = alphas, npar = 'alpha', legend = T, sex = 'man')

