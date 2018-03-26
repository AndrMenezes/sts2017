rm(list = ls())
setwd('C:/Users/User/Dropbox/4° Série/Tópicos Especiais em Estatística/Trabalhos/Trabalho 4/Scripts' )
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}
printmrow <- function(x) cat(cat(x,sep=" & "),"\\\\ \n")

filesemvs  <- list.files(pattern = 'emvs')[c(2, 1)]
filesbayes <- list.files(pattern = 'bayes')[c(2, 1)]


for(i in 1:length(filesemvs))
{
  if(i == 1) aux <- c(4, 8)
  if(i == 2) aux <- c(3, 6)
  emvs         <- read.delim(filesemvs[i])[-aux, ]
  bayes        <- read.delim(filesbayes[i])
  emvs$emv     <- paste0(FF(emvs$Estimate), " (", FF(emvs$StandardError), ")")
  bayes$mean   <- paste0(FF(bayes$Mean), " (", FF(bayes$StdDev), ")")
  emvs$ic      <- paste0("(", FF(emvs$Lower), ", " ,FF(emvs$Upper), ")")
  bayes$hpd    <- paste0("(", FF(bayes$HPDLower), ", ", FF(bayes$HPDUpper), ")")
  k            <- length(unique(emvs$Parameter))
  
  mat           <- cbind(emvs[, -c(3, 4, 5, 6)], bayes[, -c(1, 2, 3, 4, 5, 6)])
  mat$Parameter <- paste0("$\\", mat$Parameter, "$")
  mat$gender    <- c(paste0("\\multirow{", k, "}{*}{", mat$gender[1] , "}"), rep("", k -1), 
                     paste0("\\multirow{", k, "}{*}{", mat$gender[k+1] , "}"), rep("", k - 1))
  mat           <- as.matrix(mat[, c(1, 2, 3, 5, 4, 6)])
  sink(file = 'resumo-app2.tex', append = T)
  if(i == 1)
  {
    cat('\\begin{table}[H]', "\n")
    cat(paste0("\\caption{Estimativas de máxima verossimilhança e resumos posterioris.}"), "\n")
    cat(paste0("\\centering \n \\begin{tabular}{lccccc} \\toprule \n Gênero & Parâmetro & EMV (E.P.) & Média (D.P.) & IC $95\\%$ & HPD $95\\%$ \\\\ \n"))
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