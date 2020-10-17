## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----table-simple, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
df <- cbind(all = c("$\\textbf{all}$","$\\text{AvgRelA}^{[m]}_j$","$\\vdots$","$\\text{AvgRelA}^{[k]}_j$","$\\vdots$",
                         "$\\text{AvgRelA}^{[1]}_j$","$\\vdots$","$\\text{AvgRelA}_j$"),
                 uts = c("$\\textbf{uts}$","$\\text{AvgRelA}^{[m]}_{a,j}$","$\\vdots$","$\\text{AvgRelA}^{[k]}_{a,j}$","$\\vdots$",
                         "$\\text{AvgRelA}^{[1]}_{a,j}$","$\\vdots$","$\\text{AvgRelA}_{a,j}$"),
                 bts = c("$\\textbf{bts}$","$\\text{AvgRelA}^{[m]}_{b,j}$","$\\vdots$","$\\text{AvgRelA}^{[k]}_{b,j}$","$\\vdots$",
                         "$\\text{AvgRelA}^{[1]}_{b,j}$","$\\vdots$","$\\text{AvgRelA}_{b,j}$"))
rownames(df) <- c("","$\\textbf{m}$","$\\vdots$","$\\textbf{k}$","$\\vdots$","$\\textbf{1}$","$\\vdots$","$\\textbf{all}$")
knitr::kable(df,align='cccc',escape = F, col.names = rep("",3))

## ----table-simple2, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
df <- cbind(col1 = c("$\\textbf{m}$","$\\textbf{1}$","$\\text{RelA}^{[m],1}_{1,j}$","$\\vdots$",
                         "$\\text{RelA}^{[m],1}_{i,j}$","$\\vdots$","$\\text{RelA}^{[m],1}_{n,j}$"),
            col2 = c("$\\dots$","$\\dots$","$\\dots$","","$\\dots$","","$\\dots$"),
            col3 = c("", "$\\textbf{1}$","$\\text{RelA}^{[k],1}_{1,j}$","$\\vdots$",
                         "$\\text{RelA}^{[k],1}_{i,j}$","$\\vdots$","$\\text{RelA}^{[k],1}_{n,j}$"),
            col4 = c("$\\textbf{k}$","$\\dots$","$\\dots$","","$\\dots$","","$\\dots$"),
            col5 = c("", "$\\mathbf{h_k}$","$\\text{RelA}^{[k],h_k}_{1,j}$","$\\vdots$",
                         "$\\text{RelA}^{[k],h_k}_{i,j}$","$\\vdots$","$\\text{RelA}^{[k],h_k}_{n,j}$"),
            col6 = c("$\\dots$","$\\dots$","$\\dots$","","$\\dots$","","$\\dots$"),
            col7 = c("", "$\\textbf{1}$","$\\text{RelA}^{[1],1}_{1,j}$","$\\vdots$",
                         "$\\text{RelA}^{[1],1}_{i,j}$","$\\vdots$","$\\text{RelA}^{[1],1}_{n,j}$"),
            col8 = c("$\\textbf{1}$","$\\dots$","$\\dots$","","$\\dots$","","$\\dots$"),
            col9 = c("", "$\\mathbf{m}$","$\\text{RelA}^{[1],m}_{1,j}$","$\\vdots$",
                         "$\\text{RelA}^{[1],m}_{i,j}$","$\\vdots$","$\\text{RelA}^{[1],m}_{n,j}$"))
rownames(df) <- c("${\\cal K}$","$\\textbf{h}$","$\\textbf{1}$","$\\vdots$","$\\textbf{i}$","$\\vdots$","$\\textbf{n}$")
knitr::kable(df,align='ccccccccc',escape = F, col.names = rep("",9))

## ----table-simple3, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
df <- cbind(col1 = c("$\\textbf{m}$","$\\text{AvgRelA}^{[m]}_{1,j}$","$\\vdots$",
                         "$\\text{AvgRelA}^{[m]}_{i,j}$","$\\vdots$","$\\text{AvgRelA}^{[m]}_{n,j}$"),
            col2 = c("$\\dots$","$\\dots$","","$\\dots$","","$\\dots$"),
            col3 = c("$\\textbf{k}$", "$\\text{AvgRelA}^{[k]}_{1,j}$","$\\vdots$",
                         "$\\text{AvgRelA}^{[k]}_{i,j}$","$\\vdots$","$\\text{AvgRelA}^{[k]}_{n,j}$"),
            col4 = c("$\\dots$","$\\dots$","","$\\dots$","","$\\dots$"),
            col5 = c("$\\textbf{1}$","$\\text{AvgRelA}^{[1]}_{1,j}$","$\\vdots$",
                         "$\\text{AvgRelA}^{[1]}_{i,j}$","$\\vdots$","$\\text{AvgRelA}^{[1]}_{n,j}$"),
            col6 = c("$\\dots$","$\\dots$","","$\\dots$","","$\\dots$"),
            col7 = c("$\\mathbf{all}$","$\\text{AvgRelA}_{1,j}$","$\\vdots$",
                         "$\\text{AvgRelA}_{i,j}$","$\\vdots$","$\\text{AvgRelA}_{n,j}$"))
rownames(df) <- c("","$\\textbf{1}$","$\\vdots$","$\\textbf{i}$","$\\vdots$","$\\textbf{n}$")
knitr::kable(df,align='ccccccc',escape = F, col.names = rep("",7))

## ----table-simple4, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
df <- cbind(col1 = c("$\\textbf{m}$","$\\textbf{1}$","$\\text{AvgRelA}^{[m],1}_{j}$",
                         "$\\text{AvgRelA}^{[m],1}_{a,j}$","$\\text{AvgRelA}^{[m],1}_{b,j}$"),
            col2 = c("$\\dots$","$\\dots$","$\\dots$","$\\dots$","$\\dots$"),
            col3 = c("", "$\\textbf{1}$","$\\text{AvgRelA}^{[k],1}_{j}$",
                         "$\\text{AvgRelA}^{[k],1}_{a,j}$","$\\text{AvgRelA}^{[k],1}_{b,j}$"),
            col4 = c("$\\textbf{k}$","$\\dots$","$\\dots$","$\\dots$","$\\dots$"),
            col5 = c("", "$\\mathbf{h_k}$","$\\text{AvgRelA}^{[k],h_k}_{j}$",
                         "$\\text{AvgRelA}^{[k],h_k}_{a,j}$","$\\text{AvgRelA}^{[k],h_k}_{b,j}$"),
            col6 = c("$\\dots$","$\\dots$","$\\dots$","$\\dots$","$\\dots$"),
            col7 = c("", "$\\textbf{1}$","$\\text{AvgRelA}^{[1],1}_{j}$",
                         "$\\text{AvgRelA}^{[1],1}_{a,j}$","$\\text{AvgRelA}^{[1],1}_{b,j}$"),
            col8 = c("$\\textbf{1}$","$\\dots$","$\\dots$","$\\dots$","$\\dots$"),
            col9 = c("", "$\\mathbf{m}$","$\\text{AvgRelA}^{[1],m}_{j}$",
                         "$\\text{AvgRelA}^{[1],m}_{a,j}$","$\\text{AvgRelA}^{[1],m}_{b,j}$"))
rownames(df) <- c("${\\cal K}$","$\\textbf{h}$","$\\textbf{all}$","$\\textbf{a}$","$\\textbf{b}$")
knitr::kable(df,align='ccccccccc',escape = F, col.names = rep("",9))

## ----eval=FALSE---------------------------------------------------------------
#  library(FoReco)
#  data(FoReco_data)
#  
#  # Cross-temporal framework
#  oct_recf <- octrec(FoReco_data$base, m = 12, C = FoReco_data$C,
#                     comb = "bdshr", res = FoReco_data$res)$recf
#  oct_score <- score_index(recf = oct_recf,
#                           base = FoReco_data$base,
#                           test = FoReco_data$test, m = 12, nb = 5)
#  
#  # Cross-sectional framework#'
#  # monthly base forecasts
#  id <- which(simplify2array(strsplit(colnames(FoReco_data$base), split = "_"))[1, ] == "k1")
#  mbase <- t(FoReco_data$base[, id])
#  # monthly test set
#  mtest <- t(FoReco_data$test[, id])
#  # monthly residuals
#  id <- which(simplify2array(strsplit(colnames(FoReco_data$res), split = "_"))[1, ] == "k1")
#  mres <- t(FoReco_data$res[, id])
#  # monthly reconciled forecasts
#  mrecf <- htsrec(mbase, C = FoReco_data$C, comb = "shr", res = mres)$recf
#  # score
#  hts_score <- score_index(recf = mrecf, base = mbase, test = mtest, nb = 5)
#  
#  # Temporal framework
#  data(FoReco_data)
#  # top ts base forecasts ([lowest_freq' ...  highest_freq']')
#  topbase <- FoReco_data$base[1, ]
#  # top ts residuals ([lowest_freq' ...  highest_freq']')
#  topres <- FoReco_data$res[1, ]
#  # top ts test ([lowest_freq' ...  highest_freq']')
#  toptest <- FoReco_data$test[1, ]
#  # top ts recf ([lowest_freq' ...  highest_freq']')
#  toprecf <- thfrec(topbase, m = 12, comb = "acov", res = topres)$recf
#  # score
#  thf_score <- score_index(recf = toprecf, base = topbase, test = toptest, m = 12)

