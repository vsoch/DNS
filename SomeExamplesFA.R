#####################################
# Some examples for factor analysis
# Prepared by Esther Salazar
# Oct. 17
#####################################

source("functions_quest_analysis.R")

library(qgraph)
library(psych)
library(pheatmap)

# ***********************************
# NEO data and groups 
# ***********************************
data(big5)
data(big5groups)

# ***********************************
# Factor model 
# ***********************************
# Number of factors
K <- 5
big5efa <-factanal(big5, factors=K, rotation="promax", scores="regression")

# Extracting factor loadings 
big5loadings <- loadings(big5efa)

# Plot for factor loadings (nice one)
qgraph.loadings(big5loadings,groups=big5groups,rotation="promax",minimum=0.2,cut=0.4,vsize=c(1,10),borders=FALSE,vTrans=200)

# Image for factor loadings 
pheatmap(big5loadings)

# Plot to identify specific variables for each factor
fa.diagram(big5loadings,rsize=0.5,cut=0.1,,main="factor loadings")

aux <- likefadiagram(big5loadings)
bin.load <- aux$binary
load <- aux$load; load.NA <- aux$load.NA

# Extracting especific variables for each factor
NEOAClab <- rep(NA,240)
labs <- c("Neuroticism", "Extraversion", "Openness", "Agreeableness", "Conscientiousness")
for(i in 1:5){
	ss <- seq(i,240,by=5)
	NEOAClab[ss] <- labs[i]
}

for(ix in 1:K){
	NEOcat.F <- NEOAClab[c(1:240)[which(bin.load[,ix]==1)]]
	table <- as.matrix(sort(table(NEOcat.F),decreasing=T))
	print( cbind(round(table,3) , round(table/sum(table),3)) )
	print( sum(bin.load[,ix]) )
}

