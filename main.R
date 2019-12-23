## Web: https://archive.ics.uci.edu/ml/datasets/Water+Treatment+Plant

pkg <- c("tidyverse", "class", "MASS", "VIM", "mice", "pheatmap",
         "GGally", "ggfortify", "plotly", "factoextra", "FactoMineR",
         "corrplot", "ggbiplot", "gridExtra", "vegan", "coRanking", "clValid")
lapply(pkg, library, character.only = TRUE)

#####################################
#--- Read files & Data Wrangling ---#
#####################################

VAR <- "./data/water-treatment.names" %>%
  scan(what = "", sep = "\n", skip = 64, n = 38) %>%
  trimws %>% strsplit("\\s{2, }") %>%
  reduce(rbind) %>% unname %>%
  as.data.frame(stringsAsFactors = F) %>%
  mutate(V2 = gsub("-", "_", V2))
# names with "-" will make imputing process error

df <- "./data/water-treatment.data" %>%
  read.csv(header = F, na.strings = "?", stringsAsFactors = F) %>%
  mutate(V1 = as.Date(gsub("D-", "", V1), "%d/%m/%y")) %>%
  set_names(c("Date", VAR$V2)) %>%
  arrange(Date) %>% column_to_rownames("Date")

water <- kNN(as.matrix(df), k = 5, numFun = mean)[, 1:38]
water <- data.frame(water, row.names = row.names(df))

####################
#--- Data Label ---#
####################

lab.file <- dir("./data/label", full.names = T)

LAB <- lab.file %>% map(function(file){
  scan(file, what = "", sep = ",") %>% trimws %>% 
    strsplit("to") %>% map(function(x){
      y <- as.Date(gsub("D-", "", trimws(x)), "%d/%m/%y")
      if(length(y) == 2) {
        y <- seq(y[1], y[2], by = "1 day")
      } ; return(y)
    }) %>% reduce(c)
}) %>% set_names(c("Normal_Low_Influent", "Normal_over_Mean", "Normal",
                   "Settler_Problem", "Solid_Overload", "Storm")) %>%
  enframe %>% unnest %>% set_names(c("Label", "Date"))

water2 <- water %>% rownames_to_column("Date") %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(LAB)

water2$Label[is.na(water2$Label)] <-
  knn(water2[!is.na(water2$Label), ][-c(1, 40)],
      water2[ is.na(water2$Label), ][-c(1, 40)],
      water2$Label[!is.na(water2$Label)], k = 5) %>% as.character

########################
#--- Missing values ---#
########################

aggr(df, numbers = T, only.miss = T, sortVars = T,
     cex.lab = 1.2, cex.axis = 0.7, gap = 1,
     col = RColorBrewer::brewer.pal(9, "Set1")[c(2, 5, 8)])

na.split1 <- na.omit(df)
na.split2 <- df[rowSums(is.na(df)) > 0, ]
layout.mat <- matrix(c(rep(1, 3), rep(2, 2), 3, rep(4, 3)), 3, 3)
layout(layout.mat)
matrixplot(df, cex.axis = 0.5, interactive = F, main = "Water-Treatment")
matrixplot(na.split1, cex.axis = 0.5, interactive = F, main = "Complete Data")
matrixplot(na.split2, cex.axis = 0.5, interactive = F, main = "Missing Value")
matrixplot(water, cex.axis = 0.5, interactive = F, main = "Imputed by KNN with K=5")

##################
#--- Imputing ---#
##################

nxp.df <- prod(dim(df))
missing.prop <- sum(is.na(df)) / nxp.df
df.complete <- as.matrix(na.omit(df))
nxp.df.complete <- prod(dim(df.complete))
set.seed(100)
ij <- sample(1:nxp.df.complete, floor(nxp.df.complete * missing.prop))
df.complete.na <- df.complete
df.complete.na[ij] <- NA

SSE.knn <- list(KNN.Mean = mean, KNN.Median = median,
                KNN.Trimmean = function(x){mean(x, trim = 0.2)}) %>%
  map( function(.x) {
    SSE.knn <- numeric(0)
    for(i in 1:10){
      imp.knn <- kNN(df.complete.na, k = i, numFun = .x)[, 1:38]
      SSE.knn <- c(SSE.knn, sum((df.complete - imp.knn)^2)) }
    return(SSE.knn)
  } ) %>% as.data.frame %>% t

SSE.mice <- c(mean = "mean", pmm = "pmm", norm = "norm") %>%
  map( function(.x) {
    imp <- mice(df.complete.na, m = 5, maxit = 5, method = .x, seed = 100) %>%
      complete(action = 2)
    return(sum((df.complete - imp) ^ 2))
  } ) %>% unlist

# mean: Unconditional mean imputation
# pmm: Predictive mean matching
# norm: Bayesian linear regression

sort(c(knn = min(SSE.knn), SSE.mice))

########################################
#--- Exploratory Data Analysis(EDA) ---#
########################################

## Kernel Density
water %>% select(-Q_E, -ZN_E, -contains("COND"), -starts_with("RD")) %>%
  filter(water2$Label %in% c("Normal", "Normal_Low_Influent", "Normal_over_Mean")) %>%
  gather(Variable, Value) %>%
  mutate(Process = str_sub(Variable, -1) %>% factor(levels = c("E", "P", "D", "S")),
         Index = str_split(Variable, "_") %>% map(~ first(.)) %>% unlist) %>%
  ggplot(aes(x = Value, fill = Process)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Index, nrow = 2, scales = "free") +
  theme(text = element_text(face = 4))


## Boxplot
water %>% select(-Q_E, -ZN_E, -contains("COND"), -starts_with("RD")) %>%
  filter(water2$Label %in% c("Normal", "Normal_Low_Influent", "Normal_over_Mean")) %>% 
  gather(Variable, Value) %>%
  mutate(Process = str_sub(Variable, -1) %>% factor(levels = c("E", "P", "D", "S")),
         Index = str_split(Variable, "_") %>% map(~ first(.)) %>% unlist) %>% 
  ggplot(aes(x = Process, y = Value, colour = Index)) +
  geom_boxplot() +
  facet_wrap(~ Index, nrow = 1, scales = "free") +
  theme(legend.position = "none",
        text = element_text(face = 4))

## Correlation Coefficient Heatmap
water %>% cor %>% round(2) %>% as.table() %>% as.data.frame() %>%
  ggplot(aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Freq), colour = "grey50") +
  geom_text(aes(label = Freq), size = 2.5, fontface = 2,
            data = . %>% filter(abs(Freq) >= 0.4)) +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_distiller(palette = "Spectral")

water %>% cor %>% round(2) %>%
  pheatmap(color = colorRampPalette(rev(fields::larry.colors()))(100),
           display_numbers = matrix(ifelse(abs(.) >= 0.4, ., ""), nrow(.)),
           fontsize_number = 7.5, treeheight_col = 0)

## Parallel Coordinate Plot
water3 <- rbind(water2[water2$Label == "Normal", ][sample(1:330, 200), ],
                water2[water2$Label != "Normal", ]) # sample the normal observations

p1 <- water3[-1] %>% filter(Label %in% c("Normal", "Settler_Problem", "Solid_Overload", "Storm")) %>%
  select(contains("_E"), Label) %>%
  ggparcoord(1:(ncol(.)-1), ncol(.)) +
  ggtitle("Process1 : Input to Plant") + theme_bw() +
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.position = c(0.11, 0.82),
        legend.text = element_text(size = 8, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "gray50"))

p2 <- water3[-1] %>% filter(Label %in% c("Normal", "Settler_Problem", "Solid_Overload", "Storm")) %>%
  select(contains("_P"), -starts_with("RD"), Label) %>%
  ggparcoord(1:(ncol(.)-1), ncol(.)) +
  ggtitle("Process2 : Input to Primary Settler") + theme_bw() +
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.position = "none")

p3 <- water3[-1] %>% filter(Label %in% c("Normal", "Settler_Problem", "Solid_Overload", "Storm")) %>%
  select(contains("_D"), -starts_with("RD"), Label) %>%
  ggparcoord(1:(ncol(.)-1), ncol(.)) +
  ggtitle("Process3 : Input to Secondary Settler") + theme_bw() +
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.position = "none")

p4 <- water3[-1] %>% filter(Label %in% c("Normal", "Settler_Problem", "Solid_Overload", "Storm")) %>%
  select(contains("_S"), -starts_with("RD"), Label) %>%
  ggparcoord(1:(ncol(.)-1), ncol(.)) +
  ggtitle("Process4 : Output") + theme_bw() +
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.position = "none")

grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)

##################################################
#------------ n > p Dataset : 527x38 ------------#
##################################################

scale.water <- scale(water) 
scale.water <- ifelse(scale.water >= 5, 5, scale.water)
scale.water <- ifelse(scale.water <= -5, -5, scale.water)

#####################
#--- n > p : PCA ---#
#####################

water.class <- water2[, 40]
class.names <- unique(water.class)
unique.class <- as.integer(as.factor(class.names))

water.cov.pca <- PCA(scale.water, scale.unit = FALSE, graph = FALSE)
fviz_eig(water.cov.pca, addlabels = TRUE, ylim = c(0, 30))

cov.pca.var <- get_pca_var(water.cov.pca)
var.kms <- kmeans(cov.pca.var$contrib, centers = 3, nstart = 25)
kms.grp <- as.factor(var.kms$cluster)
fviz_pca_var(water.cov.pca, col.var = kms.grp, palette = c("blue", "green", "red"), legend.title = "Cluster")
fviz_pca_var(water.cov.pca, col.var = "contrib", gradient.cols = c("black", "yellow", "red4"))

cov.pca.dim1 <- water.cov.pca$ind$coord[, 1]
cov.pca.dim2 <- water.cov.pca$ind$coord[, 2]
cov.pca.dim3 <- water.cov.pca$ind$coord[, 3]
col <- c("black", "red", "blue", "green", "pruple", "gray")[unique.class]
cov.pca.plot <- plot_ly(as.data.frame(scale.water), x = ~ cov.pca.dim1,
                        y = ~ cov.pca.dim2, z = ~ cov.pca.dim3, 
                        color = ~as.factor(water2[, 40]) ,
                        colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = '1st PCA Component'),
                      yaxis = list(title = '2nd PCA Component'),
                      zaxis = list(title = '3nd PCA Component')))
cov.pca.plot
water.cov.pca <- PCA(scale.water, scale.unit = FALSE, graph = FALSE)
HCPC(water.cov.pca, nb.clust = 4)

########################
#--- n > p : Isomap ---#
########################

par(mfrow = c(2, 2))
k.range <- c(3, 5, 7, 9)
for (i in 1:4){
  water.isomap <- isomap(dist(scale.water), ndim = 2, k = k.range[i])
  plot(water.isomap, col = as.factor(water2[,40]), main = paste("ISOMAP with k =", k.range[i]))
  legend("topright", class.names, pch = 1, col = unique.class, ncol = 1, x.intersp = 0.3)
}

#######################
#--- LCMC Function ---#
#######################

cm.UL_K <- function(K, n) {
  tmp <- matrix(F, n, n)
  tmp[1:K, 1:K] <- T
  tmp
}

myLCMC <- function (Q, K = 1:nrow(Q)) {
  nQ <- nrow(Q)
  nK <- length(K)
  N <- nQ + 1
  if (nK < 0.2 * nQ) {
    lcmc <- numeric(nK)
    for (i in 1:nK) {
      k <- K[i]
      lcmc[i] <- k/(1 - N) + sum(Q[cm.UL_K(k, nQ)])/N/k
    }
  }
  else {
    lcmc_ <- diag(apply(apply(Q, 2, cumsum), 1, cumsum))/(1:nQ)/N - (1:nQ)/nQ
    lcmc <- lcmc_[K]
  }
  lcmc
}

############################
#--- n > p : LCMC - PCA ---#
############################

pca <- princomp(scale.water)
Q.pca <- coranking(scale.water, pca$score[,1:2], input = "data")
imageplot(Q.pca)
lcmc.pca <- LCMC(Q.pca, K = 5:10)

water.isomap <- isomap(dist(scale.water), ndim=2, k=9)
Q.iso <- coranking(scale.water, water.isomap$points[,1:2], input = "data")
lcmc.isomap <- LCMC(Q.iso, K = 5:10)

lcmc.pca.index <- c()
for (i in 2:7){
  Q.pca <- coranking(scale.water, pca$score[,1:i], input = "data")
  lcmc.pca <- LCMC(Q.pca, K = 3:10)
  lcmc.pca.index <- rbind(lcmc.pca.index,lcmc.pca)
}

row.names(lcmc.pca.index) <- c(paste("PCA with",2:7,"Comp"))
colnames(lcmc.pca.index) <- c(paste("K =", 3:10))
lcmc.pca.index

###############################
#--- n > p : LCMC - ISOMAP ---#
###############################

lcmc.iso.index <- c()
for (i in 2:7){
  water.isomap <- isomap(dist(scale(water)), ndim=i, k=4)
  Q.iso <- coranking(scale.water, water.isomap$points[,1:2], input = "data")
  lcmc.isomap <- LCMC(Q.iso, K = 5:10)
  lcmc.iso.index <- rbind(lcmc.iso.index,lcmc.isomap)
}
row.names(lcmc.iso.index) <- c(paste("ISOMAP with K =",2:7))
colnames(lcmc.iso.index) <- c(paste("K =", 5:10))
lcmc.iso.index

############################
#--- n > p : Clustering ---#
############################

cl.internal <- clValid(scale(water), 2:6, clMethods=c("hierarchical"), 
                       validation=c("internal"), method = c("complete"))
summary(cl.internal)
cl.stability <- clValid(scale(water), 2:6, clMethods=c("hierarchical"), 
                        validation=c("stability"), method = "complete")
summary(cl.stability)

ann_row <- data.frame(Label = water2$Label, row.names = row.names(water))
hm.c <- pheatmap(scale(water), annotation_row = ann_row,
                 cluster_cols = F, show_rownames = F, show_colnames = T,
                 clustering_method = "complete", cutree_rows = 2,
                 main = "Complete Linkage with Heatmap of Water Treatment data")
hm5.c <- pheatmap(scale.water, annotation_row = ann_row,
                  cluster_cols = F, show_rownames = F, show_colnames = T,
                  clustering_method = "complete", cutree_rows = 2,
                  main = "Complete Linkage with upper/lower bound of Water Treatment")
gridExtra::grid.arrange(hm.c[[4]], hm5.c[[4]], nrow = 1)

hm.s <- pheatmap(scale(water), annotation_row = ann_row,
                 cluster_cols = F, show_rownames = F, show_colnames = T,
                 clustering_method = "single", cutree_rows = 3,
                 main = "Single Linkage with Heatmap of Water Treatment data")
hm5.s <- pheatmap(scale.water, annotation_row = ann_row,
                  cluster_cols = F, show_rownames = F, show_colnames = T,
                  clustering_method = "single", cutree_rows = 8,
                  main = "Single Linkage with upper/lower bound of Water Treatment")
gridExtra::grid.arrange(hm.s[[4]], hm5.s[[4]], nrow = 1)

hm.a <- pheatmap(scale(water), annotation_row = ann_row,
                 cluster_cols = F, show_rownames = F, show_colnames = T,
                 clustering_method = "average", cutree_rows = 2,
                 main = "Average Linkage with Heatmap of Water Treatment data")
hm5.a <- pheatmap(scale.water, annotation_row = ann_row,
                  cluster_cols = F, show_rownames = F, show_colnames = T,
                  clustering_method = "average", cutree_rows = 2,
                  main = "Average Linkage with upper/lower bound of Water Treatment")
gridExtra::grid.arrange(hm.a[[4]], hm5.a[[4]], nrow = 1)

#################################################
#------------ n = p Dataset : 38x38 ------------#
#################################################

n <- dim(water)[1]
p <- dim(water)[2]
label.pro <- which(water2[, 40] == "Settler_Problem" |
                   water2[, 40] == "Solid_Overload" |
                   water2[,40] == "Storm")

set.seed(1234)
selected.id <- sample(1:n, 24)
id <- sort(c(label.pro, selected.id))
water.sub <- water[id, ]
dim(water.sub)

#####################
#--- n = p : PCA ---#
#####################

mycol <- as.factor(water2[id,40])
autoplot(prcomp(water.sub, scale = T), data = water2[id, ], colour = 'Label', frame = T)
fviz_eig(PCA(water.sub, scale = T), addlabels = TRUE, ylim = c(0, 30), main = "Variances of Principal Components")
fviz_contrib(PCA(water.sub, scale = T), choice = "var", axes = 1:2, top=10)
fviz_pca_var(PCA(water.sub, scale = T), col.var = "contrib",
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE)

########################
#--- n = p : Isomap ---#
########################

mycol.2 <- mycol
levels(mycol.2) <- c('coral1','goldenrod3','green3','turquoise2','royalblue1','mediumorchid1')

par(mfrow=c(1,3))
water.isomap <- isomap(dist(scale(water.sub)), ndim=2, k=3) # try different k
plot(water.isomap, col=as.vector(mycol.2), pch=16,xlab="Isomap.dim1", ylab="Isomap.dim2",main="Isomap with K=3")
legend("bottomright", legend=levels(mycol),col=levels(mycol.2), pch=16)

water.isomap <- isomap(dist(scale(water.sub)), ndim=2, k=5) # try different k
plot(water.isomap, col=as.vector(mycol.2), pch=16,xlab="Isomap.dim1", ylab="Isomap.dim2",main="Isomap with K=5")
legend("bottomright", legend=levels(mycol),col=levels(mycol.2), pch=16)

water.isomap <- isomap(dist(scale(water.sub)), ndim=2, k=7) # try different k
plot(water.isomap, col=as.vector(mycol.2), pch=16,xlab="Isomap.dim1", ylab="Isomap.dim2",main="Isomap with K=7")
legend("bottomright", legend=levels(mycol),col=levels(mycol.2), pch=16)

############################
#--- n = p : LCMC - PCA ---#
############################

pca <- princomp(scale(water.sub))
d <- matrix(nrow=6, ncol=6)
row.names(d) <- as.character(2:7)
colnames(d) <- paste("K=",5:10)
for (i in 2:7){
  Q.pca <- coranking(scale(water.sub), pca$scores[,1:i], input = "data")
  lcmc.pca <- LCMC(Q.pca, K = 5:10)
  lcmc.pca
  d[i-1,] <- lcmc.pca 
}
d

###############################
#--- n = p : LCMC - ISOMAP ---#
###############################

d2 <- matrix(nrow=6, ncol=6)
row.names(d2) <- as.character(2:7)
colnames(d2) <- paste("K=",5:10)
for (i in 2:7){
  water.isomap <- isomap(dist(scale(water.sub)), ndim=i,k=3)
  Q.pca <- coranking(scale(water.sub), water.isomap$points, input = "data")
  lcmc.pca <- LCMC(Q.pca, K = 5:10)
  lcmc.pca
  d2[i-1,] <- lcmc.pca
}
d2

############################
#--- n = p : Clustering ---#
############################

water.col <- as.data.frame(water2[id,40])
rownames(water.col) <- rownames(water.sub)
colnames(water.col) <- "Label"
pheatmap(scale(water.sub), clustering_method = "complete", annotation_row = water.col)

####################################
#--- n = p : Cluster Validation ---#
####################################

cv.i <- clValid(scale(water.sub), 2:6, clMethods=c("hierarchical"), validation="internal", method = "complete") 
cv.s <- clValid(scale(water.sub), 2:6, clMethods=c("hierarchical"), validation="stability", method = "complete")
summary(cv.i)
summary(cv.s)

pheatmap(scale(water.sub), clustering_method = "complete", annotation_row = water.col, cutree_rows = 6, main="Complete Linkage")
pheatmap(scale(water.sub), clustering_method = "single", annotation_row = water.col, cutree_rows = 2, main="Single Linkage")
pheatmap(scale(water.sub), clustering_method = "average", annotation_row = water.col, cutree_rows = 3, main="Average Linkage")

#################################################
#------------ n < p Dataset : 20x38 ------------#
#################################################

set.seed(12345)
ind <- grep("^Normal", water2$Label)
selected.id1 <- sample(ind, 13)
selected.id2 <- sample((1:527)[-ind], 7)
water20 <- water[c(selected.id1, selected.id2), ]
group <- water2$Label[c(selected.id1, selected.id2)]

############################
#--- n < p : PCA & LCMC ---#
############################

pca <- prcomp(water20, scale. = T)

fviz_eig(pca, addlabels = TRUE, ylim = c(0, 35), main = "")

ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = group,
         ellipse = TRUE, circle = TRUE) +
  theme(legend.position = c(0.9, 0.9),
        legend.title = element_blank())

PCA.LCMC <- map(2:7, function(x){
  Q.pca <- coranking(water20, pca$x[, 1:x], input = "data")
  lcmc.pca <- LCMC(Q.pca, K = 5:10)
  return(lcmc.pca)
}) %>% reduce(rbind) %>% data.frame %>%
  set_names(paste0("K = ", 5:10)) %>% remove_rownames %>%
  mutate(Dim = 2:7) %>% dplyr::select(Dim, everything())

###############################
#--- n < p : Isomap & LCMC ---#
###############################

ISOMAP.LCMC <- map(2:7, function(x){
  iso <- isomap(dist(water20), ndim = x, k = 3)
  Q.iso <- coranking(water20, iso$points, input = "data")
  lcmc.iso <- LCMC(Q.iso, K = 5:10)
  return(lcmc.iso)
}) %>% reduce(rbind) %>% data.frame %>%
  set_names(paste0("K = ", 5:10)) %>% remove_rownames %>%
  mutate(Dim = 2:7) %>% dplyr::select(Dim, everything())

par(mfrow = c(2, 3))
seq(3, 13, 2) %>% map(function(x){
  iso <- isomap(dist(water20), ndim = 2, k = x)
  plot(iso, col = as.integer(factor(group)) + 1,
       pch = 16, main = paste0("Isomap with k = ", x))
  legend("top", legend = unique(group), pch = 16, col = 2:5, ncol = 2)
})

####################################
#--- n < p : Cluster Validation ---#
####################################

int.complete <- clValid(water20, nClust = 3:8, clMethods = "hierarchical",
                        method = "complete", validation = "internal")
stab.complete <- clValid(water20, nClust = 3:8, clMethods = "hierarchical",
                         method = "complete", validation = "stability")
summary(int.complete)
summary(stab.complete)

#-------

int.single <- clValid(water20, nClust = 3:8, clMethods = "hierarchical",
                      method = "single", validation = "internal")
stab.single <- clValid(water20, nClust = 3:8, clMethods = "hierarchical",
                       method = "single", validation = "stability")
summary(int.single)
summary(stab.single)

#-------

int.average <- clValid(water20, nClust = 3:8, clMethods = "hierarchical",
                       method = "average", validation = "internal")
stab.average <- clValid(water20, nClust = 3:8, clMethods = "hierarchical",
                        method = "average", validation = "stability")
summary(int.average)
summary(stab.average)

#-------

int.ward <- clValid(water20, nClust = 3:8, clMethods = "hierarchical",
                    method = "ward", validation = "internal")
stab.ward <- clValid(water20, nClust = 3:8, clMethods = "hierarchical",
                     method = "ward", validation = "stability")
summary(int.ward)
summary(stab.ward)

############################
#--- n < p : Clustering ---#
############################

water20.scale <- scale(water20)
ann_row <- data.frame(Label = group, row.names = row.names(water20))
hm1 <- pheatmap(water20.scale, annotation_row = ann_row,
                show_rownames = F, show_colnames = F, treeheight_col = 0,
                clustering_method = "complete", cutree_rows = 3,
                main = "Complete Linkage")
hm2 <- pheatmap(water20.scale, annotation_row = ann_row,
                show_rownames = F, treeheight_col = 0,
                clustering_method = "single", cutree_rows = 3,
                main = "Single Linkage")
hm3 <- pheatmap(water20.scale, annotation_row = ann_row,
                show_rownames = F, show_colnames = F, treeheight_col = 0,
                clustering_method = "average", cutree_rows = 7,
                main = "Average Linkage")
hm4 <- pheatmap(water20.scale, annotation_row = ann_row,
                show_rownames = F, treeheight_col = 0,
                clustering_method = "ward.D", cutree_rows = 3,
                main = "Ward's Method")

gridExtra::grid.arrange(hm1[[4]], hm2[[4]], nrow = 2)
gridExtra::grid.arrange(hm3[[4]], hm4[[4]], nrow = 2)
