library(MASS)
library(tidyverse)
library(RColorBrewer)
library(gplots)

setwd("C:/Users/Rafael Ruggiero/Documents/Pos Doc/Paper Doc/e-Life/Dados Abertos OSF/Codes")

df <- read.csv2(file ='Mv_matrix.csv',sep = , header=TRUE, na.strings = '..')
df$Grupo <- as.factor(df$Grupo)  




df2 <- mutate_all(df[,7:ncol(df)], function(x) as.numeric(as.character(x)))
ind <- sapply(df2, is.numeric)

#Scaling the numerical values
df2[ind] <- lapply(df2[ind], scale)

######################################################################
#Heatmap

cc = c('#99cccc','#99cccc','#99cccc','#99cccc','#99cccc','#99cccc',"#E495A5","#E495A5","#E495A5","#E495A5","#E495A5","#E495A5","#E495A5","#E495A5","#E495A5")

#convert to matrix
mtscaled <- as.matrix(df2)


# Dissimilarity matrix
d <- dist(df2, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "ward.D" )


# transpose the matrix and cluster columns
hc.cols <- hclust(dist(t(mtscaled)))

## draw heatmap 
colMain <- colorRampPalette(brewer.pal(8, "RdBu"))(25)
heatmap(mtscaled, Rowv=as.dendrogram(hc1), Colv=as.dendrogram(hc.cols), RowSideColors=cc, scale='none', col=rev(colMain))


