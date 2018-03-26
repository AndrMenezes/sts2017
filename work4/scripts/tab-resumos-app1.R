rm(list = ls())
setwd('C:/Users/User/Dropbox/4° Série/Tópicos Especiais em Estatística/Trabalhos/Trabalho 4/Scripts' )
source('C:/Users/User/Dropbox/4° Série/Tópicos Especiais em Estatística/Trabalhos/Trabalho 4/Scripts/MOEW.R' )
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}
printmrow <- function(x) cat(cat(x,sep=" & "),"\\\\ \n")

filesemvs  <- list.files(pattern = 'emvs')[-c(1, 2)]
filesbayes <- list.files(pattern = 'bayes')[-c(1, 2)]

n <- c(100, 300, 500)
for(i in 1:length(filesbayes))
{
  emvs         <- read.delim(filesemvs[i])
  bayes        <- read.delim(filesbayes[i])
  emvs$emv     <- paste0(FF(emvs$Estimate), " (", FF(emvs$StandardError), ")")
  bayes$mean   <- paste0(FF(bayes$Mean), " (", FF(bayes$StdDev), ")")
  emvs$ic      <- paste0("(", FF(emvs$Lower), ", " ,FF(emvs$Upper), ")")
  bayes$hpd    <- paste0("(", FF(bayes$HPDLower), ", ", FF(bayes$HPDUpper), ")")
  
  mat           <- cbind(emvs[, -c(2, 3, 4, 5)], bayes[, -c(1, 2, 3, 4, 5)])
  mat$Parameter <- paste0("$\\", mat$Parameter, "$")
  mat           <- as.matrix(mat[, c(1, 2, 4, 3, 5)])
  aux           <- c(paste0("\\multirow{3}{*}{", n[i], "}"), "", "")
  mat           <- cbind(aux, mat)

  sink(file = 'resumo-app1.tex', append = T)
  if(i == 1)
  {
    cat('\\begin{table}[H]', "\n")
    cat(paste0("\\caption{Estimativas de máxima verossimilhança e resumos posterioris.}"), "\n")
    cat(paste0("\\centering \n \\begin{tabular}{lccccc} \\toprule \n $n$ & Parâmetro & EMV (E.P.) & Média (D.P.) & IC $95\\%$ & HPD $95\\%$ \\\\"))
  }
  invisible(apply(mat, 1, printmrow))
  if(i != length(filesbayes))  cat('\\hdashline ')
  if(i == length(filesbayes))
  {
    cat('\\bottomrule \n')
    cat('\\end{tabular}', "\n") 
    cat('\\end{table}')
  }
  sink()
}

################################## Histograma ##################################
dados <- read.delim('sim-MOEW.txt')[, -c(1, 5)]
head(dados)
x1 <- dados[which(dados[,i] == 1),]$ti
x2 <- dados[which(dados$n2 == 1),]$ti
x3 <- dados[which(dados$n3 == 1),]$ti


for(i in 1:length(filesbayes))
{
  x      <- dados[which(dados[,i] == 1),]$ti
  h      <- hist(x, plot = F)
  xs     <- seq(0, max(h$breaks), l = 10000)
  emvs   <- read.delim(filesemvs[i])$Estimate
  bayes  <- read.delim(filesbayes[i])$Mean
  femv   <- dmoeweibull(x = xs, lambda = emvs[1],  beta = emvs[2],  alpha = emvs[3])
  fbayes <- dmoeweibull(x = xs, lambda = bayes[1], beta = bayes[2], alpha = bayes[3])
  Ry     <- range(femv)
  Rx     <- range(xs)
  pdf(paste0('hist-', i, '.pdf'), width = 8)
  par(mar = c(3.2, 3.2, 1.0, 1.0), cex = 1.8)
  hist(x, probability = T, col = '#BDBDBD', ylim = Ry, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', main = '')
  box();rug(x)
  lines(xs, femv,   col = "#E41A1C", lwd = 2)  
  lines(xs, fbayes, col = "#0080ff", lwd = 2)  
  axis(side = 1, at = seq(0, Rx[2], l = 5),     labels = FF(seq(0, Rx[2], l = 5), 1),     padj = -0.5)
  axis(side = 2, at = seq(Ry[1], Ry[2], l = 5), labels = FF(seq(Ry[1], Ry[2], l = 5), 1), padj = -0.5)
  mtext(text = 'x', side = 1, line = 2, cex = 1.8)
  mtext(text = 'Densidade', side = 2, line = 2.2, cex = 1.8)
  if(i == 1)
  {
    title('n = 100')
    legend('topright', legend = c('EMV', 'Bayes'), col = c("#E41A1C", "#0080ff"), lwd = 2, bty = 'n', inset = 0.03)
  }
  if(i == 2) title('n = 300')
  if(i == 3) title('n = 500')
  graphics.off()
}



library( RColorBrewer)
brewer.pal(10, name = 'Greys')
RColorBrewer::display.brewer.all()

N <- M - burin
thin = 100
m <- seq(0, N, thin)
length(m)



