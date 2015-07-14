library(devtools)
source('~/oxford/switch/sigmoidr/R/switch.R')
library(gridExtra)
library(cowplot)
library(dplyr)
library(embeddr)

# let’s try on some monocle data ------------------------------------------

library(monocle)
data(HSMM)
cds <- HSMM[, HSMM$State %in% 1:2]

mrf_genes <- c('CDK1', 'MEF2C', 'MYH3', 'MYOG', 'ID1')
mrf_long_genes <- row.names(fData(HSMM))[match(mrf_genes, fData(HSMM)$gene_short_name)]

mrf_indices <- match(mrf_long_genes, featureNames(cds))
plot_genes_in_pseudotime(cds[mrf_indices, ])

mrf_expr <- exprs(cds)[mrf_indices,]
pst <- cds$Pseudotime
x <- mrf_expr[3,]

model <- fitModel(x, pst, 'nb')
null_model <- fitModel(x, 'nb', sigmoidal = FALSE)
plot_model(model, x, pst)
diffExprTest(x, t, dist = 'nb')

pvals_nb <- apply(mrf_expr, 1, diffExprTest, pst, dist = 'nb')
pvals_norm <- apply(log2(mrf_expr+1), 1, diffExprTest, pst, dist = 'normal')

plots <- apply(mrf_expr, 1, function(y) {
  y <- log2(y+1) # sqrt(e)
  alt_model <- fitModel(y, t, dist = 'normal')
  print(alt_model$par)
  plt <- plot_model(alt_model, y, pst) #+ ggtitle(mrf_genes[i]) #+ scale_y_log10()
  plt
})

plot_grid(plotlist = plots, labels = mrf_genes)


load("~/oxford/embeddr/data/sce_pst.Rdata")
fpkm(sce_23_kept) <- 10^(exprs(sce_23_kept) - 1)
X <- exprs(sce_23_kept)
pt <- pseudotime(sce_23_kept)

pvals <- apply(X, 1, diffExprTest, pt, 'normal')

df_pval <- data.frame(pval = pvals)
df_pval$gene <- row.names(df_pval)
df_pval$qval <- p.adjust(df_pval$pval, method = 'BH')

is_de <- df_pval$qval < 0.01

pst_models <- apply(fpkm(sce_23_kept[is_de,]), 1, fitModel, t = pseudotime(sce_23_kept), dist = 'nb')
pst_hypers <- fitPstModelHypers(pst_models)

null_models <- apply(fpkm(sce_23_kept[!is_de,]), 1, fitModel, dist = 'nb', sigmoidal = FALSE)
null_hypers <- fitNullModelHypers(null_models)
