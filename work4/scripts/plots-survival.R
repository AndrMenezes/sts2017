rm(list = ls())
library(survival)
source('C:/Users/User/Dropbox/4° Série/Tópicos Especiais em Estatística/Competing Risk/Ajustes/Turnbull.R')
wd1   <- 'C:/Users/User/Dropbox/4° Série/Tópicos Especiais em Estatística/Trabalhos/Trabalho 4/Scripts' 
wd2   <- 'C:/Users/User/Dropbox/4° Série/Tópicos Especiais em Estatística/Trabalhos/Trabalho 4' 
FF    <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}
smoew <- function(x, mu, beta, alpha) {alpha * exp(-(x / mu)^beta) / (1 - (1 - alpha) * exp(-(x / mu)^beta))}

setwd(wd1)
filesemvs  <- list.files(pattern = 'emvs')[c(2, 1)]
filesbayes <- list.files(pattern = 'bayes')[c(2, 1)]
dados      <- read.delim('diabetes.txt')

fit     <- survfit(Surv(time = left, time2 = right, type = 'interval2') ~ gender, data = dados)
t1      <- seq(0, 44, l = 1000)
t2      <- seq(0, 40, l = 1000)

bayes1 <- read.delim(filesbayes[1])
emvs1  <- read.delim(filesemvs[1])[-c(4, 8),]
bayes2 <- read.delim(filesbayes[2])
emvs2  <- read.delim(filesemvs[2])[-c(3, 6),]

moew.mle.fem   <- smoew(x = t1, mu = emvs1$Estimate[1], beta = emvs1$Estimate[2], alpha = emvs1$Estimate[3])
moew.bayes.fem <- smoew(x = t1, mu = bayes1$Mean[1],    beta = bayes1$Mean[2],    alpha = bayes1$Mean[3])
moew.mle.mas   <- smoew(x = t2, mu = emvs1$Estimate[4], beta = emvs1$Estimate[5], alpha = emvs1$Estimate[6])
moew.bayes.mas <- smoew(x = t2, mu = bayes1$Mean[4],    beta = bayes1$Mean[5],    alpha = bayes1$Mean[6])

my_plot1 <- function(gender)
{
  if(gender == 'fem')
  {
    t      <- t1
    km     <- fit[1, ]
    Semvs  <- moew.mle.fem
    Sbayes <- moew.bayes.fem
    nome   <- 'Feminino'
  }
  if(gender == 'man')
  {
    t      <- t2
    km     <- fit[2, ]
    Semvs  <- moew.mle.mas
    Sbayes <- moew.bayes.mas
    nome   <- 'Masculino'
  }
  R <- range(t)
  pdf(file = paste0('km-', gender, '.pdf'), width = 8)
  par(mar = c(3.2, 3.2, 1.0, 1.0), cex = 1.8)
  plot(km, mark.time = TRUE, conf.int = FALSE, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', lwd = 2, lty = 2)
  axis(side = 1, at = seq(R[1], R[2], l = 5), FF(seq(R[1], R[2], l = 5), 0), hadj = 0.65, padj = -0.4)
  axis(side = 2, at = seq(0, 1, l = 5), FF(seq(0, 1, l = 5), 2), padj = 0.4)
  mtext(text = 'Tempo',                  side = 1, line = 2, cex = 1.8)
  mtext(text = 'Sobrevivência estimada', side = 2, line = 2, cex = 1.8)
  lines(t, Semvs,  col = 2, lwd = 2)
  lines(t, Sbayes, col = 4, lwd = 2)
  legend('topright', legend = c('Turnbull', 'EMV', 'Bayes'), col = c(1, 2, 4), lwd = 2, lty = c(2, 1, 1), inset = 0.02, bty = 'n')
  title(main = nome, cex = 1.8)
  graphics.off()
}

setwd(wd2)
my_plot1(gender = 'fem')
my_plot1(gender = 'man')


my_plot2 <- function(gender)
{
  if(gender == 'fem')
  {
    t      <- t1
    km     <- fit[1, ]
    Smoew  <- moew.mle.fem
    Sweib  <- weib.mle.fem
    nome   <- 'Feminino'
  }
  if(gender == 'man')
  {
    t      <- t2
    km     <- fit[2, ]
    Smoew  <- moew.mle.mas
    Sweib  <- weib.mle.mas
    nome   <- 'Masculino'
  }
  R <- range(t)
  pdf(file = paste0('km-', gender,  '-2.pdf'), width = 8)
  par(mar = c(3.2, 3.2, 1.0, 1.0), cex = 1.8)
  plot(km, mark.time = TRUE, conf.int = FALSE, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', lwd = 2, lty = 2)
  axis(side = 1, at = seq(R[1], R[2], l = 5), FF(seq(R[1], R[2], l = 5), 0), hadj = 0.65, padj = -0.4)
  axis(side = 2, at = seq(0, 1, l = 5), FF(seq(0, 1, l = 5), 2), padj = 0.4)
  mtext(text = 'Tempo',                  side = 1, line = 2, cex = 1.8)
  mtext(text = 'Sobrevivência estimada', side = 2, line = 2, cex = 1.8)
  lines(t, Smoew,  col = 2, lwd = 2)
  lines(t, Sweib,  col = 4, lwd = 2)
  legend('topright', legend = c('Turnbull', 'MOEW', 'Weibull'), col = c(1, 2, 4), lwd = 2, lty = c(2, 1, 1), inset = 0.02, bty = 'n')
  title(main = nome, cex = 1.8)
  graphics.off()
}

my_plot2(gender = 'fem')
my_plot2(gender = 'man')






# When the data set includes left censored or interval censored data (or both),
# then the EM approach of Turnbull is used to compute the overall curve. When the baseline method
# is the Kaplan-Meier, this is known to converge to the maximum likelihood estimate.



tobs <- fit[1, ]$time
pemp <- fit[1, ]$surv
moew.mle.fem   <- smoew(x = tobs, mu = emvs$Estimate[1], beta = emvs$Estimate[2], alpha = emvs$Estimate[3])
moew.bayes.fem <- smoew(x = tobs, mu = bayes$Mean[1],    beta = bayes$Mean[2],    alpha = bayes$Mean[3])
mean( (pemp - moew.bayes.fem)^2 )
mean( (pemp - moew.mle.fem)^2 )

moew.mle.mas   <- smoew(x = t2, mu = emvs$Estimate[4], beta = emvs$Estimate[5], alpha = emvs$Estimate[6])
moew.bayes.mas <- smoew(x = t2, mu = bayes$Mean[4],    beta = bayes$Mean[5],    alpha = bayes$Mean[6])
