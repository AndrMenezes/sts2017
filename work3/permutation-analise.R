rm(list = ls(all.names = TRUE))
setwd('C:/Users/User/Dropbox/4° Série/Tópicos Especiais em Estatística/Trabalhos/Trabalho 3')
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}
library(labestData)
library(ggplot2)
library(reshape2)
library(nlme)
library(lme4)
library(lmPerm)
library(goftest)
library(hnp)

# Hotelling's T2 Test -----------------------------------------------------

T2.Hotelling <-function(X1, X2, M)
{
  p  <- ncol(X1)
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  n  <- n1 + n2
  T2 <- function(X1, X2, perm = FALSE)
  {
    if(perm == TRUE)
    {
      Z  <- rbind(X1, X2)
      Z  <- apply(Z, 2, sample, replace = FALSE)
      X1 <- Z[1:n1, ]
      X2 <- Z[(n1+1):n, ]
    }
    M1   <- apply(X1, 2, mean)
    S1   <- var(X1)
    M2   <- apply(X2, 2, mean)
    S2   <- var(X2)
    S    <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n - 2) 
    stat <- ((n1 * n2) / n) * t(M1 - M2) %*% solve(S) %*% (M1 - M2)
    return(as.numeric(stat))
  }
  T2.orig <- T2(X1, X2)
  F0.orig <- (n - p - 1) / (p * (n -2)) * T2.orig
  T2.perm <- replicate(M, T2(X1, X2, perm = TRUE))
  F0.perm <- (n - p - 1) / (p * (n -2)) * T2.perm
  p.asym  <- 1 - pf(F0.orig, p, n - p - 1)
  p.perm  <- (1 + sum(T2.perm >= T2.orig)) / (M + 1)
  output  <- list(statistic = c(T2 = T2.orig, F0 = F0.orig), p.value = c(asymptotic = p.asym, permutation = p.perm), 
                  T2.perm = T2.perm, F0.perm = F0.perm)
  return(output)
}


# Nearest neighbor test ---------------------------------------------------

NN.test <- function(X, Y, k=NULL, M)
{
  Z    <- rbind(X, Y)
  n1   <- nrow(X)
  n    <- nrow(Z)
  if(is.null(k)) k <- n - 1
  Tn.k <- function(Z, k)
  {
    NN   <- FNN::get.knn(Z, k = k)$nn.index
    TN.k <- (sum(NN[1:n1, ] < n1 + 0.5) + sum(NN[(n1 + 1):n, ] > n1 + 0.5)) / (n * k)
    return(TN.k)
  }
  TNk.orig <- Tn.k(Z, k)
  TNk.perm <- replicate(n = M, expr = Tn.k(Z = apply(Z, 2, sample, replace = FALSE), k = k))
  p.value  <- (1 + sum(TNk.perm >= TNk.orig)) / (M + 1)
  output   <- list(statistic = TNk.orig, p.value = p.value, TNk.perm = TNk.perm) 
  return(output)
}

# Examples ----------------------------------------------------------------

## Example on simulated variables
set.seed(1212)
X1 <- mvtnorm::rmvnorm(n = 10, mean = rep(0, 4), sigma = diag(4))
X2 <- mvtnorm::rmvnorm(n = 10, mean = rep(5, 4), sigma = diag(4))

outT2 <- T2.Hotelling(X1, X2, M = 1000000)
outNN <- NN.test(X1, X2, k = 6, M = 1000000)

set.seed(1212)
X1 <- rmvt(n = 10, sigma = diag(4), df = 10)
X2 <- rmvt(n = 10, sigma = diag(4), df = 10)

outT2 <- T2.Hotelling(X1, X2, M = 1000)
outNN <- NN.test(X1, X2, k = 6, M = 1000)


## Manly Table 1.1 - page 4
X1    <- subset(ManlyTb1.1, ManlyTb1.1$sobrev == 'S')[, -6]
X2    <- subset(ManlyTb1.1, ManlyTb1.1$sobrev == 'N')[, -6]
outT2 <- T2.Hotelling(X1, X2, M = 10000)
outNN <- NN.test(X1, X2, k = 10, M = 10000)

pdf('T2perm_1.pdf', width = 9)
par(mar = c(3.0, 3.0, 1.0, 1.0), cex = 1.8)
hist(outT2$T2.perm, probability = TRUE, col = '#0080ff', breaks = 20, border = 'orange', ylim = c(0, 0.14),
     xlim = range(outT2$T2.perm), xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', main = ''); box()
axis(side = 1, at = seq(min(outT2$T2.perm), max(outT2$T2.perm), l = 5), FF(seq(min(outT2$T2.perm), max(outT2$T2.perm), l = 5), 2))
axis(side = 2, at = seq(0, 0.14, l = 5), FF(seq(0, 0.14, l = 5), 2))
mtext(text = 'Densidade', side = 2, line = 2.2, cex = 1.8)
abline(v = outT2$statistic[1], col = 2, lwd = 2)
legend('topright', inset = 0.01, bty = 'n', col = c(2, NA), lwd = c(2, NA), 
       legend = c(parse(text = paste0('T^2 ==', FF(outT2$statistic[1], 4))), 
                  parse(text = paste0('hat(p) ==', FF(outT2$p.value[2], 4)))))
graphics.off()

pdf('TNkperm_1.pdf', width = 9)
par(mar = c(3.0, 3.0, 1.0, 1.0), cex = 1.8)
hist(outNN$TNk.perm, probability = TRUE, col = '#0080ff', breaks = 20, border = 'orange', ylim = c(0, 17),
     xlim = range(outNN$TNk.perm), xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', main = ''); box()
axis(side = 1, at = seq(min(outNN$TNk.perm), max(outNN$TNk.perm), l = 5), FF(seq(min(outNN$TNk.perm), max(outNN$TNk.perm), l = 5), 2))
axis(side = 2, at = seq(0, 17, l = 5), FF(seq(0, 17, l = 5), 2))
mtext(text = 'Densidade', side = 2, line = 2.2, cex = 1.8)
abline(v = outNN$statistic, col = 2, lwd = 2)
legend('topright', inset = 0.01, bty = 'n', col = c(2, NA), lwd = c(2, NA), 
       legend = c(parse(text = paste0('T[10] ==', FF(outNN$statistic, 4))), 
                  parse(text = paste0('hat(p) ==', FF(outNN$p.value, 4)))))
graphics.off()

## Package Hotelling
data(container.df)
split.data <- split(container.df[,-1],container.df$gp)
X1         <- split.data[[1]]
X2         <- split.data[[2]]
outT2      <- T2.Hotelling(X1, X2, M = 10000)
outNN      <- NN.test(X1, X2, k = 6, M = 10000)


## Book Hand and Taylor
mat   <- read.table('salsolinol2.txt', sep = ',', header = T)
X1    <- mat[mat$grupo == 1,][, -1]
X2    <- mat[mat$grupo == 2,][, -1]
outT2 <- T2.Hotelling(X1, X2, M = 100000)
outNN <- NN.test(X1, X2, k = 10, M = 100000)
outNN$p.value

var(X1)
var(X2)
apply(X1, 2, mean)
apply(X2, 2, mean)

pdf('T2perm_2.pdf', width = 9)
par(mar = c(3.0, 3.0, 1.0, 1.0), cex = 1.8)
hist(outT2$T2.perm, probability = TRUE, col = '#0080ff', breaks = 20, border = 'orange', ylim = c(0, 0.085),
     xlim = range(outT2$T2.perm), xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', main = ''); box()
axis(side = 1, at = seq(min(outT2$T2.perm), max(outT2$T2.perm), l = 5), FF(seq(min(outT2$T2.perm), max(outT2$T2.perm), l = 5), 2))
axis(side = 2, at = seq(0, 0.085, l = 5), FF(seq(0, 0.08, l = 5), 2))
mtext(text = 'Densidade', side = 2, line = 2.2, cex = 1.8)
abline(v = outT2$statistic[1], col = 2, lwd = 2)
legend('topright', inset = 0.01, bty = 'n', col = c(2, NA), lwd = c(2, NA), 
       legend = c(parse(text = paste0('T^2 ==', FF(outT2$statistic[1], 4))), 
                  parse(text = paste0('hat(p) ==', FF(outT2$p.value[2], 4)))))
graphics.off()


ymax <- max(hist(outNN$TNk.perm, probability = TRUE)$density)
pdf('TNkperm_2.pdf', width = 9)
par(mar = c(3.0, 3.0, 1.0, 1.0), cex = 1.8)
hist(outNN$TNk.perm, probability = TRUE, col = '#0080ff', breaks = 20, border = 'orange', ylim = c(0, ymax),
     xlim = range(outNN$TNk.perm), xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', main = ''); box()
axis(side = 1, at = seq(min(outNN$TNk.perm), max(outNN$TNk.perm), l = 5), FF(seq(min(outNN$TNk.perm), max(outNN$TNk.perm), l = 5), 2))
axis(side = 2, at = seq(0, ymax, l = 5), FF(seq(0, ymax, l = 5), 2))
mtext(text = 'Densidade', side = 2, line = 2.2, cex = 1.8)
abline(v = outNN$statistic, col = 2, lwd = 2)
legend('topright', inset = 0.01, bty = 'n', col = c(2, NA), lwd = c(2, NA), 
       legend = c(parse(text = paste0('T[8] ==', FF(outNN$statistic, 4))), 
                  parse(text = paste0('hat(p) ==', FF(outNN$p.value, 4)))))
graphics.off()


mat       <- cbind(ind = 1:14, mat)
df        <- melt(mat, id.vars = c('ind', 'grupo'))
names(df) <- c('ind', 'grupo', 'dia', 'salsolinol')
df$grupo  <- as.factor(df$grupo)
df$dia    <- rep(1:4, each = 14) 

ggplot(data = df, aes(x = dia, y = salsolinol, col = grupo, group = grupo)) +
  stat_summary(fun.y = mean, geom="point", size = 2, shape = 15) +
  stat_summary(fun.y = mean, geom="line") 

ggplot(df, aes(x=dia, y=salsolinol, col = grupo)) +
  geom_point() + facet_wrap(~ind) + geom_smooth(method = 'lm', se = F) + theme(legend.position="top")


levels(df$grupo) <- c('Moderada', 'Severa')
df$dia.f         <- as.factor(df$dia)
levels(df$dia.f) <- c('1° dia', '2° dia', '3° dia', '4° dia')

ggplot(df, aes(y = salsolinol, x = grupo, fill = grupo)) +
  facet_wrap(~dia.f) +
  geom_boxplot() +
  labs(y = 'Salsolinol (mmol)', x = '', fill = 'Dependência') +
  theme(axis.title.x    = element_blank(),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(), 
        text            = element_text(size = 14), 
        legend.position = 'top') -> pp
ggsave(filename = 'bp2.pdf', plot = pp, device = 'pdf', width = 9, height = 7)


reg <- aov(salsolinol ~ dia.f + grupo + dia.f*grupo, data = df)
anova(reg)
res.std <- rstandard(reg)
shapiro.test(res.std)
ks.test(res.std, 'pnorm')
cvm.test(res.std, 'pnorm')
ad.test(res.std, 'pnorm')


hnp.plot <- function(myhnp, nome)
{
  Rx <- range(myhnp$x)
  Ry <- c(min(myhnp$lower), max(myhnp$residuals))
  pdf(file = paste0("hnp-", nome, ".pdf"), width = 11, height = 7)
  par(mar = c(3.2, 3.2, 1.5, 1.5), cex = 1.8)
  plot(myhnp$x, myhnp$residuals, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', cex = 0.3, pch = 19,
       xlim = Rx, ylim = Ry)
  lines(myhnp$x, myhnp$lower); lines(myhnp$x, myhnp$upper); lines(myhnp$x, myhnp$median)
  mtext("Percentil da N(0, 1)", side = 1, line = 2.0, cex = 1.8)
  mtext("Resíduos", side = 2, line =2, cex = 1.8)
  abline(h=seq(Ry[1], Ry[2], l = 5), v=seq(Rx[1], Rx[2], l = 5), col = "gray", lty = "dotted")
  axis(1, seq(Rx[1], Rx[2], l = 5), FF(seq(Rx[1], Rx[2], l = 5), 2))
  axis(2, seq(Ry[1], Ry[2], l = 5), FF(seq(Ry[1], Ry[2], l = 5), 2))
  graphics.off()
}
myhnp <- hnp(reg, halfnormal = F)
hnp.plot(myhnp, 'normal')


bartlett.test(salsolinol~ dia.f, data = df)
bartlett.test(salsolinol~ grupo, data = df)
df$dg <- with(df, factor(paste(dia.f, grupo, sep = '-')))
bartlett.test(salsolinol~ dg, data = df)

tapply(X = df$salsolinol, df$dia.f, mean )
xtable(anova(reg), digits = 4)

re.lm <- lmer(salsolinol ~ grupo * dia.f + (1|ind), data = df) 
summary(re.lm)
anova(re.lm)
hnp(re.lm)

M  <- 10000
F0 <- anova(reg)$`F value`[-4]
FP <- replicate(M, anova(lm(salsolinol ~ sample(dia.f) + sample(grupo) + sample(dia.f)*sample(grupo), data = df))$`F value`[-4])
(1 + rowSums(FP >= F0)) / (M + 1)

hist.perm <- function(x, lab, t0, df1, df2, curva = FALSE)
{
  M      <- length(x)
  pvalue <- (1 + sum(x >= t0)) / (M + 1)
  R      <- range(x)
  ymax   <- max(hist(x, probability = TRUE)$density)

  pdf(paste0(lab, '.pdf'), width = 9)
  par(mar = c(3.0, 3.0, 1.0, 1.0), cex = 1.8)
  hist(x, probability = TRUE, col = '#0080ff', breaks = 20, border = 'orange', ylim = c(0, ymax),
       xlim = R, xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', main = ''); box()
  axis(side = 1, at = seq(R[1], R[2], l = 5), FF(seq(R[1], R[2], l = 5), 2))
  axis(side = 2, at = seq(0, ymax, l = 5), FF(seq(0, ymax, l = 5), 2))
  mtext(text = 'Densidade', side = 2, line = 2.2, cex = 1.8)
  abline(v = t0, col = 2, lwd = 2)
  legend('topright', inset = 0.01, bty = 'n', col = c(2, NA), lwd = c(2, NA), 
         legend = c(parse(text = paste0('F[0] ==', FF(t0, 4))), 
                                parse(text = paste0('hat(p) ==', FF(pvalue, 4)))))
  if(curva == TRUE)
  {
    curve(expr = stats::df(x = x, df1 = df2, df2 = df2), from = 0, to = R[2], add = TRUE, col = 'black', lwd = 2)
  }
  graphics.off()
}

hist.perm(x = FP[1, ], lab = 'dia-perm', t0 = F0[1], df1 = 3, df2 = 48, curva = FALSE)
hist.perm(x = FP[2, ], lab = 'grupo-perm', t0 = F0[2], df1 = 1, df2 = 48, curva = FALSE)
hist.perm(x = FP[3, ], lab = 'dg-perm', t0 = F0[3], df1 = 3, df2 = 48, curva = FALSE)

summary(aovp(salsolinol ~ dia.f*grupo, data = df), perm = 'Exact')




# Comparando mais de dois grupos ------------------------------------------
df      <- melt(ManlyTb1.2, id.vars = "grup")
df$grup <- as.factor(df$grup)
levels(df$grup) <- c('12° e 13° dinastias', 'Período ptolemaico', 'Período romano', 'Pré-dinastico antigo',
                     'Pré-dinastico primitivo')
ggplot(df, aes(y = value, x = grup, fill = grup)) +
  facet_wrap(~variable, scales = 'free') +
  geom_boxplot() +
  labs(y = '', x = '', fill = '') +
  theme(axis.title.x    = element_blank(),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(), 
        text            = element_text(size = 14), 
        legend.position = 'top') -> pp
ggsave(filename = 'bp2.pdf', plot = pp, device = 'pdf', width = 9, height = 7)

library(car)
M      <- 10000
mtests <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
trat   <- as.factor(ManlyTb1.2$grup)
Y      <- as.matrix(ManlyTb1.2[, -1])
mod    <- lm(Y ~ trat) 
res    <- Anova(mod) 
x <- summary(res)
x$multivariate.tests$trat$

t.orig <- summary(manova(Y ~ trat), test = mtests[2])$stats[1, 3]
t.perm <- replicate(M, summary(manova(Y ~ sample(trat, size = nrow(Y), replace = FALSE)), test = mtests[2])$stats[1, 3])
(1 + sum(t.perm >= t.orig)) / (M + 1)

summary(manova(Y ~ trat), test = mtests[2])
summary(manova(Y ~ trat), test = mtests[3])
summary(manova(Y ~ trat), test = mtests[4])









#########################################
# Teste de permutação com dados pareados
########################################

n         = nrow(x)
nperm     = 2^n
diferenca = x[,3] - x[,2]
t.obs     = mean(diferenca)
matriz.si = as.matrix(expand.grid(rep(list(c(-1,1)), n)))
matriz.st = diferenca * matriz.si
t.perm    = apply(matriz.st, 1, mean)
mean(t.perm >= t.obs)
t.test(x[,3], x[,2], paired = T)
hist(t.perm, probability = T)
curve(dnorm(x, 0, 1), -2, 2, add= T)

setwd('C:\\Users\\User\\Dropbox\\EstatisticaBayesiana\\AulasPraticas\\Aula3Bayes')
library(readxl)
dados <- data.frame(read_xls('psi2.xls', sheet = 'dados'))

n         = nrow(dados)
nperm     = 2^n
diferenca = dados[,2] - dados[,1]
t.obs     = mean(diferenca)
matriz.si = as.matrix(expand.grid(rep(list(c(-1,1)), n)))
matriz.st = diferenca * matriz.si
t.perm    = apply(matriz.st, 1, mean)
mean(t.perm >= t.obs)
t.test(dados[,2], dados[,1], paired = T)
hist(t.perm, probability = T)
curve(dnorm(x, 0, 1), -2, 2, add= T)

