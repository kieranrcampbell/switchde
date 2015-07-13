library(devtools)
source('R/normal_switch.R')
library(gridExtra)
library(cowplot)
library(dplyr)

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

pvals <- apply(log10(mrf_expr+1), 1, diff_expr_test, pst)

plots <- apply(mrf_expr, 1, function(e) {
  y <- log10(e+1) # sqrt(e)
  alt_model <- fit_alt_model(y , pst)
  plt <- plot_model(alt_model, y, pst) #+ ggtitle(mrf_genes[i]) #+ scale_y_log10()
  plt
})

plot_grid(plotlist = plots, labels = mrf_genes)

# have a look at the residuals --------------------------------------------

y <- log10(mrf_expr[3,] + 1)
alt_model <- fit_alt_model(y, pst)
predict_expr <- function(t, params) {
  L <- params[1] ; k <- params[2] ; t_0 <- params[3] ; r <- params[4]
  L / (1 + exp(-k*(t - t_0)))
}
pe <- predict_expr(pst, alt_model$par)
resid <- y - pe
qplot(resid, geom = 'density')
qplot(t, resid)





# don’t run ---------------------------------------------------------------
library(embeddr)
load("~/delete_me.Rdata")
X <- exprs(sce_23_kept)
pt <- pseudotime(sce_23_kept)

pvals <- apply(X, 1, diff_expr_test, pt)
pvals <- data.frame(t(pvals))
pvals$gene <- row.names(pvals)
pvals$qval <- p.adjust(pvals$pval, method = 'BH')

DE <- filter(pvals, qval < 0.01)

ca <- mutate(DE, grad = k)

plt_on <- ggplot(filter(ca, grad > 0.5)) + geom_point(aes(x = t0, y = grad)) +
  scale_y_log10() + ylab('Effect size') + xlab('Activation time')

plt_off <- ggplot(filter(mutate(ca, grad = -grad), grad > .5)) +
  geom_point(aes(x = t0, y = grad)) +
  scale_y_log10() + ylab('Effect size') + xlab('Activation time')

plot_grid(plt_on, plt_off, labels = c('Genes switched on', 'Genes switched off'))

df_up <- filter(ca, grad > 0.5)
df_up <- mutate(df_up, k = log10(k))
m_up <- as.matrix(select(df_up, t0, k))
kk_up <- kmeans(m_up, centers = 2)
df_up$cl <- kk_up$cluster
up_cl <- ggplot(df_up) + geom_point(aes(x = t0, y = k, color=as.factor(cl))) +
   ylab('Effect size') + xlab('Activation time')

df_down <- filter(ca, grad < -0.5)
df_down <- mutate(df_down, k = log10(abs(k)))
m_down <- as.matrix(select(df_down, t0, k))
kk_down <- kmeans(m_down, centers = 2)
df_down$cl <- kk_down$cluster
down_cl <- ggplot(df_down) + geom_point(aes(x = t0, y = k, color=as.factor(cl))) +
  ylab('Effect size') + xlab('Activation time')


plot_grid(up_cl, down_cl, labels = c('Genes switched on', 'Genes switched off'))

## look at genes up
up_1 <- df_up$gene[kk_up$cluster == 1]
up_2 <- df_up$gene[kk_up$cluster == 2]
f <- file("~/delete_me_genes.txt")
ho <- sapply(up_2, function(s) strsplit(s, '.', fixed=T)[[1]][1])
write(ho, f)
close(f)


## look at a couple of gene families
high_on <- filter(ca, grad < 1)$gene
long_gene <- fData(HSMM[high_on,])$gene_short_name

f <- file("~/delete_me_genes.txt")
ho <- sapply(high_on, function(s) strsplit(s, '.', fixed=T)[[1]][1])
write(ho, f)
close(f)

## do some gene clustering
X <- as.matrix(select(ca, t0, grad))
km <- kmeans(X, centers=5)
mm <- Mclust(X)
plot(X, col=km$cluster)

## hierarchical
switched <- ca #filter(DE, t0 < 1, t0 > 0, abs(k/L) > 5)
predicted <- apply(switched, 1, function(row) {
  params <- as.numeric(as.vector(row[2:5]))
  calc_mu(params, pseudotime(sce_23_kept))
})

colnames(predicted) <- switched$gene
d <- as.dist(1 - cor(predicted) / 2)

hc <- hclust(d, method = 'centroid')
plot(hc)

gene_classes <- cutree(hc, 4)

df_gene <- data.frame(gene=colnames(predicted), class=gene_classes)
pe <- predicted

pe <- data.frame(scale(pe)) ## scaled pst-by-gene
pe$pseudotime <- pseudotime(sce_23_kept)
## save(pe, df_gene, file='~/pe.Rdata')

pe_melted <- melt(pe, id.vars='pseudotime', value.name='expression', variable.name='gene')
pe_melted <- inner_join(pe_melted, df_gene)

## we want to plot the mean expression rather than all of it (messy)
gp_groups <- group_by(pe_melted, class, pseudotime)
mean_expr_per_group <- dplyr::summarise(gp_groups, mean_expr = mean(expression))
pe_melted <- inner_join(pe_melted, mean_expr_per_group)
## pe_melted <- arrange(pe_melted, gene, pseudotime)

ggplot(pe_melted) + geom_line(aes(x=pseudotime, y=mean_expr), color='red') +
  facet_wrap(~ class) + stat_density2d(aes(x=pseudotime, y=expression), n=150) +
  theme_bw() + ylab('Expression') # add ncol = 1 for vertical representation

## pca
M <- select(filter(DE, qval < 0.1), k, L, t0)
pca <- prcomp(M)

x <- X[100,]
qplot(pt, x)
alt_model <- fit_alt_model(x, pt)
null_model <- fit_null_model(x)
plot_model(alt_model, x, pt)

log_lik_diff <- alt_model$value - null_model$value

