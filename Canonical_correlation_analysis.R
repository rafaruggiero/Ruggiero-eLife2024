# Canonical Correlation Analsis
#
# Script for canonical correlation analysis to investigate the relationship 
# between behavioral changes and neurobiological variables
#
# Author: Rafael Naime Ruggiero 2022

library(lme4)
library(CCA) #facilitates canonical correlation analysis
library(CCP) #facilitates checking the significance of the canonical variates
library(dplyr)
library(ggpubr)

setwd("C:/Users/Rafael Ruggiero/Documents/Pos doc/Paper Doc/e-Life/Dados Abertos OSF/Codes")

# Read the CSV file into a data frame
df <- read.csv2(file ='Imputed_data_matrix.csv',sep = , header=TRUE, na.strings = '..')
df <- mutate_all(df, function(x) as.numeric(as.character(x)))

# Invert the sign of PC1
df$PC1 <- df$PC1*-1

# Separate the  behavioral and neurobiological variables
X <- df[,1:2]
Y <- df[,3:ncol(df)]

# Perform canonical correlation analysis
cc_results <- cancor(X,Y)


str(cc_results)
cc_results$xcoef
cc_results$ycoef
cc_results$cor


can_cor2 <- comput(X,Y,cc_results)
can_cor2

#test of canonical dimensions
rho <-cc_results$cor

#defining the number of observations, no of variables in first set and number of variables in second set
n <- dim(X)[1]
p <- length(X)
q <- length(Y)



#Calculating the F approximations using Wilk's Statistics
p.asym(rho, n, p, q, tstat="Wilks")
p.asym(rho, n, p, q, tstat="Hotelling")
p.perm(X, Y, nboot = 999, rhostart = 1)


# Create a data frame for plotting
dfcca <- X
dfcca$xscores1 <- can_cor2$xscores[,1]
dfcca$xscores2 <- can_cor2$xscores[,2]
dfcca$yscores1 <- can_cor2$yscores[,1]
dfcca$yscores2 <- can_cor2$yscores[,2]
dfcca$xscores1 <- can_cor2$xscores[,1]
dfcca$Grupo <- "Ctrl"
dfcca$Grupo[8:13] <- "SE"

ccascatter <- ggplot(dfcca, aes(x=xscores1*-1,y=yscores1*-1, color=Grupo, fill=Grupo))+
  geom_point(size=6) + theme_classic2() +xlab("CCX") +ylab("CCY") +  
  scale_color_manual(values=c("blue", "red"))+
  geom_abline()

ccascatter 

# barchart with added parameters
bar1 <- barplot(append(can_cor2$corr.X.xscores[,1]*-1,can_cor2$corr.Y.xscores[,1]*-1),
                main = "CCX",
                ylab = "Loadings",
                xlab = "Covariate",
                names.arg = c("PC1", "PC2","LTP","GFAP","NeuN","PV","mGLUR5"),
                col = "darkgreen",
                horiz = FALSE,
                ylim = c(-1,1))



# barchart with added parameters


bar2 <- barplot(append(can_cor2$corr.X.yscores[,1]*-1,can_cor2$corr.Y.yscores[,1]*-1),
                main = "CCY",
                ylab = "Loadings",
                xlab = "Covariate",
                names.arg = c("PC1", "PC2","LTP","GFAP","NeuN","PV","mGLUR5"),
                col = "darkorchid3",
                horiz = FALSE,
                ylim = c(-1,1))

