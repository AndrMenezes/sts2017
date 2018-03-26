rm(list = ls())
setwd('C:/Users/User/Dropbox/4° Série/Tópicos Especiais em Estatística/Regras de associação/artigo')
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}
pkg <- c('imputeTS', 'arules', 'arulesCBA', 'arulesViz' ,'class', 'MASS', 'e1071', 'vcd', 'ggplot2', 'xtable', 
         'caret', 'plyr', 'klaR')
sapply(pkg, require, character.only = T)

fl <- list.files(pattern = '.txt')

dados <- read.table(fl[2], sep = ';', skip = 16, header = T, stringsAsFactors = F)

# Formatando dados --------------------------------------------------------
names(dados) <- tolower(names(dados))
dados$x      <- NULL     
dados$data   <- as.Date(dados$data, format = '%d/%m/%Y')
ind          <- match(c("estacao", "data", "hora", "precipitacao", "tempminima"), names(dados))
zero         <- subset(dados, dados$hora == 0, select = -ind[-c(1:2)])
doze         <- subset(dados, dados$hora == 1200, select = ind[-3])
dados        <- merge(zero, doze, by = c('estacao', 'data'))
dados$mes    <- factor(format(dados$data, "%m"))
dados        <- dados[, c('data', 'tempmaxima', 'umidade.relativa.media', 'velocidade.do.vento.media', 'precipitacao', 'tempminima')]
names(dados) <- c('data', 'tempmaxima', 'umidade.relativa', 'velocidade.do.vento', 'precipitacao', 'tempminima')
head(dados)

apply(dados, 2, function(x) sum(is.na(x)))

# Decritiva ---------------------------------------------------------------
descritiva <- function(x)
{
  na <- sum(is.na(x))
  FF(c(na, min(x, na.rm = T), mean(x, na.rm = T), median(x, na.rm = T), max(x, na.rm = T), sd(x, na.rm = T),
       skewness(x, na.rm = T), kurtosis(x, na.rm = T)), 2)
}
K <- do.call(rbind, lapply(dados[, -1], descritiva))
print(xtable(K), include.rownames = T, include.colnames = F)

x <- ts(dados$tempmaxima, start = c(1961, 1), end = c(2016, 31), frequency = 31)
pdf(file = 'tmax-na.pdf', width = 8)
par(mar = c(3.2, 3.2, 1.0, 1.0), cex = 1.4)
plotNA.distribution(x = x, pch = '', main = '', ylab = '', xlab = '', xlim = c(1960, 2020))
mtext(side = 2, text = 'Temperatura máxima (°C)', line = 2.2, cex = 1.4)
graphics.off()

x <- ts(na.interpolation(dados$tempmaxima), start = c(1961, 1), end = c(2016, 31), frequency = 31)
pdf(file = 'tmax.pdf', width = 8)
par(mar = c(3.2, 3.2, 1.0, 1.0), cex = 1.4)
plotNA.distribution(x = x, pch = '', main = '', ylab = '', xlab = '', xlim = c(1960, 2020))
mtext(side = 2, text = 'Temperatura máxima (°C)', line = 2.2, cex = 1.4)
graphics.off()

plotNA.distributionBar(dados$tempmaxima)
plotNA.gapsize(dados$tempmaxima)

# Interpolação para dados faltantes ---------------------------------------
dados$tempmaxima          <- na.interpolation(dados$tempmaxima)
dados$umidade.relativa    <- na.interpolation(dados$umidade.relativa)
dados$precipitacao        <- na.interpolation(dados$precipitacao)
dados$tempminima          <- na.interpolation(dados$tempminima)
dados$velocidade.do.vento <- na.interpolation(dados$velocidade.do.vento)

dados$precipitacao <- factor(ifelse(is.na(dados$precipitacao), NA, ifelse(dados$precipitacao == 0, 'Não Chuva', 'Chuva')))

getSeason <- function(dates) 
{
  ano       <- format(dates, '%Y')
  Outono    <- as.Date(paste0(ano, "-03-20"), format = "%Y-%m-%d") 
  Inverno   <- as.Date(paste0(ano, "-06-21"), format = "%Y-%m-%d") 
  Primavera <- as.Date(paste0(ano, "-09-22"), format = "%Y-%m-%d") 
  Verao     <- as.Date(paste0(ano, "-12-21"), format = "%Y-%m-%d") 

  season <- ifelse(dates >= Outono & dates < Inverno, "Outono",
                   ifelse(dates >= Inverno & dates < Primavera, "Inverno",
                          ifelse (dates >= Primavera & dates < Verao, "Primavera", "Verão")))
  return(factor(season))
}
dados$estacao <- getSeason(dates = dados$data)
head(dados)

# Discretizando as variáveis ----------------------------------------------
dados.disc <- dados
set.seed(15)
dados.disc$tempmaxima           <- discretize(dados$tempmaxima, method = 'interval')
dados.disc$tempminima           <- discretize(dados$tempminima, method = 'interval')
dados.disc$velocidade.do.vento  <- discretize(dados$velocidade.do.vento)
dados.disc$umidade.relativa     <- discretize(dados$umidade.relativa, method = 'interval')
head(dados.disc)

# Regras de associação ----------------------------------------------------
dados.trans <- as(dados.disc[, -1], 'transactions')
summary(dados.trans)

rules <- apriori(dados.trans, 
                 parameter  = list(minlen = 3, maxlen = 4, support = 0.01, confidence = 0.7),
                 appearance = list(rhs = c("precipitacao=Chuva"), default = 'lhs'))
rules <- sort(rules, by = 'confidence')
quality(rules) <- round(cbind(quality(rules), interestMeasure(rules, transactions = dados.trans, 
                                                              measure = c('conviction', 'oddsRatio'), significance = TRUE)), 3)
inspect(rules)

subsetrules <- which(colSums(is.subset(rules, rules)) > 1) # get subset rules in vector
length(subsetrules) 
rules_rm <- rules[-subsetrules] # remove subset rules. 

inspect(rules_rm)


# Tabela  -----------------------------------------------------------------
rules_df       <- as(rules_rm, 'data.frame')[-c(1, 2, 3), ]
rules_df$rules <- as.character(rules_df$rules)
rules_df$rhs   <- sapply(1:nrow(rules_df), function(i) strsplit(x = rules_df$rules[i], split = '=>', fixed = T)[[1]][1])
rules_df$imp   <- '$\\Rightarrow$'
rules_df$lhs   <- sapply(1:nrow(rules_df), function(i) strsplit(x = rules_df$rules[i], split = '=>', fixed = T)[[1]][2])
rules_df       <- rules_df[, c('rhs', 'imp', 'lhs', 'support', 'confidence', 'lift', 'conviction', 'oddsRatio')]
rules_df       <- as.matrix(rules_df)
head(rules_df)


printmrow <- function(x) cat(cat(x,sep=" & "),"\\\\ \n")

sink(file = 'rules.tex', append = T)
cat('\\begin{landscape} \n \\begin{table}[H]', "\n")
cat(paste0("\\caption{Regras de associação.}"), "\n")
cat(paste0("\\scalefont{0.75} \\centering \n \\begin{tabular}{lcc|ccccc} \\toprule \n \\multicolumn{3}{c|}{Regra  ($ A \\Rightarrow B$)} & Suporte & Confiança & Lift & Conviction & OR \\\\ \\midrule \n"))
invisible(apply(rules_df, 1, printmrow))
cat('\\bottomrule \n')
cat('\\end{tabular}', "\n") 
cat('\\end{table} \n \\end{landscape}', "\n")
sink()


# Gráfico -----------------------------------------------------------------
x11()
plot(rules_rm, measure=c("support", "lift"), shading="confidence", cex = 1, pch = 19, xlab = '')

rules_df <- as(rules_rm, 'data.frame')
ggplot(data = rules_df, aes(x = support, y = lift, col = confidence )) +
  geom_point(size = 3.0) + labs(x = 'Suporte', y = 'Lift') +
  scale_color_gradient2(low ="white", high ="firebrick2", mid = "dodgerblue", midpoint = median(rules_df$confidence),
                       limits = c(0.7, 1.0), name = "Confiança\n",
                       breaks = c(0.7, 0.8, 0.9, 0.99)) +
  scale_y_continuous(limits = range(rules_df$lift), breaks = seq(min(rules_df$lift), max(rules_df$lift), l = 5), 
                     labels = FF(seq(min(rules_df$lift), max(rules_df$lift), l = 5), 2)) +
  scale_x_continuous(limits = range(rules_df$support), breaks = seq(min(rules_df$support), max(rules_df$support), l = 5), 
                     labels = FF(seq(min(rules_df$support), max(rules_df$support), l = 5), 2)) +
  theme(panel.grid.minor = element_blank(), text = element_text(size = 14),
        axis.title = element_text(size = 14), panel.grid.major = element_line(size = 0.5), 
        legend.position="right",legend.key.size = unit(1.0, "cm")) 
ggsave(filename = 'plotar.pdf', device = 'pdf', width = 10, height = 6)

x11()
plot(rules_rm, method = 'graph')

x11()
plot(rules_rm, method = 'grouped')

# Classificação -----------------------------------------------------------
# folds <- createFolds(y = dados$precipitacao, k = 6);str(folds)
set.seed(1212)
folds <- createDataPartition(y = dados$precipitacao, times = 999, p = 0.75)
str(folds)

# CBA ---------------------------------------------------------------------
cv_CBA <- ldply(lapply(folds, function(x)
{
  dados.treino <- as(dados.disc[x, -1], 'transactions')
  dados.teste  <- dados.disc[-x, -1]
  rules        <- apriori(dados.treino, parameter = list(minlen = 3, maxlen = 4, support = 0.01, confidence = 0.75),
                          appearance = list(rhs = c("precipitacao=Chuva","precipitacao=Não Chuva"), default = 'lhs'), 
                          control = list(verbose = F))
  rules_rm   <- rules[-which(colSums(is.subset(rules, rules)) > 1)]
  classifier <- CBA_ruleset(precipitacao ~ ., rules = rules_rm)
  pred.CBA   <- predict(classifier, dados.teste)
  df         <- data.frame(obs = dados.teste$precipitacao, pred = pred.CBA)
  med        <- c(defaultSummary(data = df, lev = levels(df$obs)), sensitivity(df$pred, df$obs))
  return(med)
}))

# KNN ---------------------------------------------------------------------
dados.s <- cbind(as.data.frame(lapply(dados[, -c(1, 5, 7)], scale)), precipitacao = dados$precipitacao, 
                 estacaoano = as.integer(dados[, 7]))
## Determinando valor de k, avaliando a acurácia
ks <- 2:12

deterk <- ldply(lapply(1:length(ks), function(i)
{
  set.seed(1212)
  amostras     <- sample(1:nrow(dados.s), size = ceiling(0.8 * nrow(dados.s)), replace = F)
  dados.treino <- dados.s[amostras, ]
  dados.teste  <- dados.s[-amostras, ]
  pred.knn     <- knn(train = dados.treino[, -5], test = dados.teste[, -5], cl = dados.treino$precipitacao, 
                      k = ks[i])
  df  <- data.frame(obs = dados.teste$precipitacao, pred = pred.knn)
  med <- c(defaultSummary(data = df, lev = levels(df$obs)), sensitivity(df$pred, df$obs), specificity(df$pred, df$obs))
  return(med)
}))

medidas <- c('Acurácia', 'Kappa', 'Sensibilidade', 'Especificidade')

for(j in 1:ncol(deterk))
{
  y      <- deterk[, j]
  maximo <- ks[which.max(y)]
  pdf(file = paste0(medidas[j], '.pdf'), width = 9)
  par(mar = c(3.2, 3.2, 0.8, 0.8), cex = 1.6)  
  plot(ks, y, type = 'o', pch = 19, xaxt = 'n', xlab = '', ylab = '', yaxt = 'n', cex = 0.9)
  points(x = maximo, y = max(y), pch = 19, cex = 0.9, col = 'red')
  R <- range(y)
  axis(side = 1, ks)
  axis(side = 2, seq(R[1], R[2], l = 5), FF(seq(R[1], R[2], l = 5), 2))
  mtext(medidas[j], side = 2, line = 2.1, cex = 1.6)
  mtext('k', side = 1, line = 2.1, cex = 1.6)
  abline(h = seq(R[1], R[2], l = 5), v = ks, col = "lightgray", lty = "dotted")
  graphics.off()
}

## Cross-Validation com k = 10
cv_KNN <- ldply(lapply(folds, function(x)
{
  dados.treino <- dados.s[x, ]
  dados.teste  <- dados.s[-x, ]
  pred.knn     <- knn(train = dados.treino[, -5], test = dados.teste[, -5], cl = dados.treino$precipitacao, 
                      k = 10)
  df  <- data.frame(obs = dados.teste$precipitacao, pred = pred.knn)
  med <- c(defaultSummary(data = df, lev = levels(df$obs)), sensitivity(df$pred, df$obs))
  return(med)
}))

# LDA ---------------------------------------------------------------------
cv_LDA <- ldply(lapply(folds, function(x)
{
  dados.treino <- dados[x, -c(1, 7)]
  dados.teste  <- dados[-x, -c(1, 7)]
  prior        <- as.numeric(prop.table(table(dados.treino$precipitacao)))
  classifier   <- lda(precipitacao ~., data = dados.treino, prior = prior)
  pred.lda     <- predict(classifier, dados.teste)$class
  df           <- data.frame(obs = dados.teste$precipitacao, pred = pred.lda)
  med <- c(defaultSummary(data = df, lev = levels(df$obs)), sensitivity(df$pred, df$obs))
  return(med)
}))


# Naive-Bayes -------------------------------------------------------------
cv_NB <- ldply(lapply(folds, function(x)
{
  dados.treino <- dados[x, -c(1, 7)]
  dados.teste  <- dados[-x, -c(1, 7)]
  # classifier <- naiveBayes(precipitacao ~ ., data = dados.treino)
  classifier   <- NaiveBayes(precipitacao ~ ., data = dados.treino, usekernel = T, kernel = 'gaussian')
  pred.nb      <- predict(classifier, dados.teste)$class
  df           <- data.frame(obs = dados.teste$precipitacao, pred = pred.nb)
  med          <- c(defaultSummary(data = df, lev = levels(df$obs)), sensitivity(df$pred, df$obs))
  return(med)
}))

resumo <- function(x)  c(rbind(apply(X = x, 2, quantile, probs = 0.025), colMeans(x), apply(X = x, 2, quantile, probs = 0.975)))
H <- rbind(resumo(x = cv_CBA[, -1]), resumo(x = cv_KNN[, -1]),
           resumo(x = cv_LDA[, -1]), resumo(x = cv_NB[, -1]))
M <- H[, c(2, 5, 8)]
H <- H[, -c(2, 5, 8)]
R <- paste0(rowSums(sapply(1:ncol(M), function(i) rank(-M[,i]))), "$^", rank(R), "$")
M <- sapply(1:ncol(M), function(i)  paste0(FF(M[, i]), "$^", rank(-M[, i], ties.method = 'min'), "$"))
H <- cbind(FF(H, 4), M)
H <- H[, c(1, 7, 2, 3, 8, 4, 5, 9, 6)]
H <- cbind(c('CBA', 'KNN', 'LDA', 'NB'), H, R)

sink(file = 'classifier.tex', append = T)
cat('\\begin{table}[H]', "\n")
cat(paste0("\\caption{Performance dos métodos de classificação.}"), "\n")
cat(paste0("\\centering \n \\begin{tabular}{lcccccccccc} \\toprule \n Método & \\multicolumn{3}{c|}{Acurácia} & \\multicolumn{3}{c|}{Kappa} & \\multicolumn{3}{c|}{Sensibilidade} & Rank\\\\ \\midrule \n"))
invisible(apply(H, 1, printmrow))
cat('\\bottomrule', '\n')
cat('\\end{tabular}', '\n') 
cat('\\end{table}')
sink()

H <- scan() 
0.7374190 0.7478762 0.7578475 0.2902622 0.3279797 0.3678279 0.9366667 0.9525407 0.9659259
0.7750374 0.7857152 0.7957150 0.4781888 0.5023196 0.5247423 0.8492593 0.8618407 0.8748148
0.7767813 0.7876304 0.7989537 0.4733772 0.4992557 0.5256240 0.8651852 0.8771672 0.8892593
0.7760214 0.7868486 0.7977080 0.4734838 0.4996328 0.5245761 0.8585000 0.8725244 0.8855556
H <- data.frame(matrix(H, nrow = 4, ncol = 9, byrow = T))  
pontual <- H[, c(2, 5, 8)]
names(pontual) <- c('Acurácia', 'Kappa', 'Sensibilidade')
met <- c('CBA', 'KNN', 'LDA', 'NB')
df  <- cbind(met, gather(pontual, 'medida', 'point'))
aux <- H[, -c(2, 5, 8)]
li  <- gather(aux, li, lower, c(1, 3, 5))[, 5]
ls  <- gather(aux, ls, upper, -c(1, 3, 5))[, 5]
df  <- cbind(df, li, ls)
ggplot(df, aes(x=medida, y = point, col = met, ymin = li, ymax = ls)) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.2) +
  labs(x = 'Medida', y = 'Valor', col = 'Método') +
  scale_y_continuous(limits = c(min(df$li), max(df$ls)), breaks = seq(min(df$li), max(df$ls), l = 5), 
                     labels = FF(seq(min(df$li), max(df$ls), l = 5), 2)) +
  theme(panel.grid.minor = element_blank(), text = element_text(size = 14),
        axis.title = element_text(size = 14), panel.grid.major = element_line(size = 0.5), 
        legend.position="right",legend.key.size = unit(1.0, "cm")) 
ggsave(filename = 'performance.pdf', width = 10, height = 6)


# Análise Descritiva ------------------------------------------------------
tabela <- function(X)
{
  nm <- names(X)
  for(i in 1:length(nm))
  {
    if(i == 1)
    {
      sink(file = 'descritiva.tex', append = T)
      cat('\\begin{table}[H]', "\n")
      cat(paste0("\\caption{Estatísticas descritivas.}"), "\n")
      cat(paste0("\\centering \n \\begin{tabular}{ccc} \\toprule \n Variável & Categorias & Frequência (\\%) \\\\ \\midrule"))
    }
    tab  <- table(X[, i], useNA = 'ifany')
    ptab <- as.matrix(prop.table(tab))[, 1]
    atab <- t(rbind(names(tab), paste0(tab, ' (', FF(100 * ptab, 3), ')')))
    for(j in 1:nrow(atab))
    {
      if(j == 1) {cat(paste0('\\multirow{', nrow(atab), '}{*}{', nvar[i], '}'), atab[j, ], sep = ' & '); cat('\\\\ \n')}
      if(j != 1 & j != nrow(atab)) {cat('', atab[j, ], sep = ' & '); cat('\\\\ \n')}
      if(j == nrow(atab)) 
      {
        cat('', atab[j, ], sep = ' & ')
        if(i == length(nm)) {cat('\\\\ \\bottomrule \n')}
        else{cat('\\\\ \\hdashline \n')} 
      }
    }
    if(i == length(nm)) {cat('\\end{tabular}', "\n"); cat('\\end{table}', "\n"); sink()}
  }
}  

nvar <- c("Precipitação", "Estação do ano" ,"Umidade relativa (\\%)", "Velocidade do vento", 
          "Temperatura máxima ($^\\circ$C)", "Temperatura mínima ($^\\circ$C)")
tabela(X = dados.disc[, c(5, 7, 3, 4, 2, 6)])

df <- do.call('rbind', lapply(3:ncol(dados.disc), 
                              function(i)
                              { 
                                nome <- names(dados.disc)[i]
                                tab  <- as.data.frame(table(dados.disc[, i]))
                                tab$Var1 <- paste0(nome, '=', tab$Var1)
                                return(tab)
                              } ))
df$por <- (df$Freq / sum(df$Freq)) * 100

ggplot(df, aes(x = reorder(Var1, -Freq), y = por)) +
  geom_bar(stat = 'identity') + labs(x = '', y = 'Percentual (%)') +
  geom_text(aes(reorder(Var1, - por), por+0.5, label = Freq)) +
  theme(axis.text.x = element_text(angle=70, hjust=1)) +
  theme(panel.grid.minor = element_blank(), text = element_text(size = 12),
        axis.title = element_text(size = 14), panel.grid.major = element_line(size = 0.5)) 
ggsave(filename = 'itemsets.pdf', device = 'pdf', width = 10, height = 6)

###### Avaliando a normalidade multivariada dos dados
X  <- dados[, -c(1, 5, 7)]
ll <- split(X, dados$precipitacao)

mh <- mahalanobis(x = ll[[1]], center = colMeans(ll[[1]]), cov = cov(ll[[1]]))

plot(density(mh, bw = 0.5), xlab = '', main = '', ylab = 'Densidade'); rug(mh)
qqplot(qchisq(ppoints(nrow(ll[[1]])), df = ncol(ll[[1]])), D2)
abline(0, 1, col = 'gray')

MVN::mardiaTest(ll[[2]])
biotools::boxM(X, dados$precipitacao)



library(caret)

set.seed(1212)
x = rnorm(20)
createFolds(x, k = 2)

folds <- createFolds(credit$default, k = 10)

cv_results <- lapply(folds, function(x) {
  credit_train <- credit[x, ]
  credit_test <- credit[-x, ]
  credit_model <- C5.0(default ~ ., data = credit_train)
  credit_pred <- predict(credit_model, credit_test)
  credit_actual <- credit_test$default
  kappa <- kappa2(data.frame(credit_actual, credit_pred))$value
  return(kappa)
})




library(ROCR)
data(ROCR.simple)
pred <- prediction(ROCR.simple$predictions, ROCR.simple$labels)
perf <- performance(pred,"tpr","fpr")
plot(perf)

library(ROCR)

help(package = 'ROCR')






