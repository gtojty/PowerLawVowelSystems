#library("ggtree")
library("ape")
library("caper")
library("phytools")
library("geiger")
library("readxl")
library("vegan")
library("motmot.2.0")
library("ggplot2")
library("phangorn")
library("phylogram")
library("diversitree")
library("picante")

setwd("G:/Structural evolution of sound system as Power Law distribution in human speech language/33 IE Languages")


################################################################################
## import phylogeny tree
tree = read.nexus("IE_temp.nex")
#plot(tree)


################################################################################
## import the trait data
data = read_xlsx("33IElanguages.xlsx", col_names = TRUE)
EffectiveDE = data$Effective_DE
FE = data$FE
x = fastBM(tree)
y = fastBM(tree)
x[1:33] = EffectiveDE[1:33]
y[1:33] = FE[1:33]
#xy = rbind(x,y)


################################################################################
## Map continuous trait evolution on the tree
## The mapping is accomplished by estimating states at internal nodes using ML
#fastAnc(tree, EffectiveDE)
fig1 = contMap(tree, x, lwd = 5, ftype = "reg", direction="rightwards")  ##Effective DE
fig2 = contMap(tree, y, lwd = 5, ftype = "reg", direction="leftwards")  ##FE



################################################################################
## Using PGLS to test the correlated evolution
## From https://www.r-phylo.org/wiki/HowTo/PGLS
# Brownian motion: the trait covariance between any pair of taxa decreases linearly with the time (in branch length) since their divergence.
Data = data.frame(log10(x), log10(y))
bm.IE = corBrownian(phy=tree)
bm.gls = gls(EffectiveDE~FE,correlation=bm.IE,data=Data)
summary(bm.gls)

# Ornstein-Uhlenbeck model: the expected covariance decreases exponentially, as governed by the parameter alpha (Martins and Hansen 1997).
ou.IE = corMartins(1,phy=tree)
ou.gls = gls(EffectiveDE~FE,correlation=ou.IE,data=Data)
summary(ou.gls)

##Alternative Choice to apply PGLS
data$log10_EffectiveDE = -log10(data$Effective_DE)
data$log10_FE = -log10(data$FE)
data = data.frame(data)
Comparative_Data = comparative.data(tree, data, IELEX_Name)
DEFE_pgls = pgls(log10_FE ~ log10_EffectiveDE, Comparative_Data, lambda='ML')
summary.pgls(DEFE_pgls)


################################################################################
## correlation and Correlated evolution
stats_cor = cor.test(EffectiveDE, FE, method = "spearman")
stats_cor

## Apply PIC to exclude the effects of IE language phylogeny, and test the correlation between two traits
pic.X = pic(-log10(x), tree, scaled = FALSE, var.contrasts=FALSE)  ##origina Effective DE
pic.Y = pic(-log10(y), tree, scaled = FALSE, var.contrasts=FALSE)  ##origina FE
stats_pic = cor.test(pic.X, pic.Y, method = "spearman")
stats_pic


################################################################################
## Calculate lambda kappa and delta for the Effective DE
x_Mat = as.matrix(x)
xkappa = transformPhylo.ML(x_Mat, phy = tree, model = "kappa",modelCIs = FALSE)
xkappa_0 = transformPhylo.ll(x_Mat, tree, model = "kappa", kappa = 0)
xkappa_1 = transformPhylo.ll(x_Mat, tree, model = "kappa", kappa = 1)

#xlambda = transformPhylo.ML(x_Mat, tree, model = "lambda", modelCIs = FALSE)
#xlambda_0 = transformPhylo.ll(x_Mat, tree, model = "lambda", lambda = 0 )
#xlambda_1 = transformPhylo.ll(x_Mat, tree, model = "lambda", lambda = 1 )

xdelta = transformPhylo.ML(x_Mat, tree, model = "delta", modelCIs = FALSE, upperBound = Inf)
xdelta_1 = transformPhylo.ll(x_Mat, tree, model = "delta", delta = 1)


## Calculate lambda kappa and delta for the FE
y_Mat = as.matrix(y)
ykappa = transformPhylo.ML(y_Mat, phy = tree, model = "kappa",modelCIs = FALSE)
ykappa_0 = transformPhylo.ll(y_Mat, tree, model = "kappa", kappa = 0)
ykappa_1 = transformPhylo.ll(y_Mat, tree, model = "kappa", kappa = 1)

#ylambda = transformPhylo.ML(y_Mat, tree, model = "lambda", modelCIs = FALSE)
#ylambda_0 = transformPhylo.ll(y_Mat, tree, model = "lambda", lambda = 0 )
#ylambda_1 = transformPhylo.ll(x_Mat, tree, model = "lambda", lambda = 1 )

ydelta = transformPhylo.ML(y_Mat, tree, model = "delta", modelCIs = FALSE, upperBound = Inf)
ydelta_1 = transformPhylo.ll(y_Mat, tree, model = "delta", delta = 1)


