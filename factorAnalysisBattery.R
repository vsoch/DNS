# FACTOR ANALYSIS ---------------------------------------------
# For each assessment, we need to do dimensionality reduction and choose
# factors that have the hightest weights to describe most of the data

# Libraries for visualizations, scree plot
library(qgraph)
library(psych)
library(pheatmap)
library(nFactors)
library(Hmisc)

# Esther Salazar functions for factor analysis
setwd("/home/vanessa/Documents/Dropbox/Code/R/DNS")
source("functions_quest_analysis.R")

# Load the battery data - note that we don't have all data for the questions defined
setwd("/home/vanessa/Documents/Work/DNS")
load("DNSBattery.Rda")

# Each group of questions is grouped based on "group" variable
uniquegroups = unique(battery$columns$GROUP)
factorAnalysis = list()
# AT END - need to deal with g=3
for (g in 1:length(uniquegroups)){
  # Vanessa - run this manually to go through questions
  group = uniquegroups[g]
  questions = battery$columns[which(battery$columns$GROUP==g),]
  data = battery$filtered[,which(colnames(battery$battery) %in% as.character(questions$LABEL))]
  
  dim(data)
  questions$LABEL
  
  # Plot the questions
  if (length(data)>0){
    par(mfrow=c(2,5))
    
    # Visualize distributions
    label = "AUDIT"
    for (d in 1:ncol(data)){
      hist(data[,d],main=colnames(data)[d],col=sample(colours(),1),xlab="")
    }
    
    # What are the scales for each question?
    cat("id","min","max","NA","tot","\n",sep="\t")
    for (d in 1:ncol(data)){
      minimum = min(data[which(!is.na(data[,d])),d])
      maximum = max(data[which(!is.na(data[,d])),d])
      numNA = length(which(is.na(data[,d])))
      total = length(data[,d])
      cat(label,minimum,maximum,numNA,total,"\n",sep="\t")
    }
      
    # Get rid of people we don't have data for
    data = data[which(as.numeric(rowSums(!is.na(data)))==ncol(data)), ]
    dim(data)
    
    # Figure out number of factors with Cattell Scree test
    # Determine Number of Factors to Extract
    ev = eigen(cor(data)) # get eigenvalues
    ap = parallel(subject=nrow(data),var=ncol(data),rep=100,cent=.05)
    nS = nScree(x=ev$values, aparallel=ap$eigen$qevpea)
    plotnScree(nS,main=paste("Scree test for",label))
    
    # Now do factor analysis
    # Do we want to normalize the weight matrix to Z scores and choose some threshold?
    K = nS$Components$noc
    cat("Number of optimal factors:",K,"\n")
    fa = factanal(data, factors=K, rotation="promax", scores="regression")  
    # Extracting factor loadings 
    loadings <- loadings(fa)
    # Visualize factor loadings before thresholding
    pheatmap(loadings,main=paste("Loadings before thresholding for",label))
    # Plot to identify specific variables for each factor
    fa.diagram(loadings,rsize=0.5,cut=0.1,,main=paste("Factor Loadings",label))
    # Binary matrix to identify the questions pointed in fa.diagram
    aux = likefadiagram(loadings)
    bin.load = aux$binary
    load = aux$load; load.NA = aux$load.NA
    # We don't want to lose our labels
    rownames(aux$load) = attr(aux$load,"dimnames")[[1]]
    colnames(aux$load) = attr(aux$load,"dimnames")[[2]]
    rownames(aux$load.NA) = attr(aux$load.NA,"dimnames")[[1]]
    colnames(aux$load.NA) = attr(aux$load.NA,"dimnames")[[2]]
    rownames(bin.load) = attr(aux$load,"dimnames")[[1]]
    colnames(bin.load) = attr(aux$load,"dimnames")[[2]]
    # Plot new heatmap
    pheatmap(bin.load,cluster_rows=T,cluster_cols=F,main=paste("Final question assignments for",label))
   
    # Save data, fa model, and factors to file
    readme = c("fa: the factor analysis of data\nF: the matrix F (N people x K factors) to be used as features\nscree: to determine number of factors\nK: number of factors\nquestions: assigned questions to each factor\ndata: raw questionnaire data with NA and missing subjects removed\nloadings: the matrix W (K factors x P questions)")
    tmp = list(fa=fa,F=fa$scores,questionnaire=label,scree=nS,K=K,questions=bin.load,data=data,loadings=loadings,README=readme)
    factorAnalysis[[label]] = tmp
    save(factorAnalysis,file="FactorAnalysisBattery.Rda")
    
    # Redo visualization stuff just for saving to PDF
    pdf(paste("distributions/",label,"_factorAnalysis.pdf",sep=""),onefile=TRUE)
    
    # Visualize factor loadings before thresholding
    pheatmap(loadings,main=paste("Loadings before thresholding for",label))
    
    # Plot new heatmap
    pheatmap(bin.load,cluster_rows=T,cluster_cols=F,main=paste("Final question assignments for",label))
    
    # The factor loadings (assigned questions)
    fa.diagram(loadings,rsize=0.5,cut=0.1,,main=paste("Factor Loadings",label))
    
    # Save the scree plot
    plotnScree(nS,main=paste("Scree test for",label))
    
    # Last, the original question scales
    for (d in 1:ncol(data)){
      hist(data[,d],main=colnames(data)[d],col=sample(colours(),1),xlab="")
    }
    dev.off()
    
    }
}

# Vanessa's original idea with Z scores
# Calculate threshold for loadings
# Z = (loadings-mean(loadings))/sd(loadings)
# hist(Z,col="orange",main=paste("Normalized Z scores for",label),xlim=c(-3,3))
# lines(x=rep(2.7,31),y=seq(0,30),col="red",lwd=2)
# lines(x=rep(-2.7,31),y=seq(0,30),col="red",lwd=2)

# Threshold the data, look at again
# Zthresh = loadings
# Zthresh[Z < 2.7] = 0
# pheatmap(Zthresh)

# Which factors will we choose?
# factors = as.numeric(which(colSums(Zthresh)>0))

# If we need to impute - FOR TESTING ONLY
# cols = unique(which(is.na(data),arr.ind=TRUE)[,2])
# for (col in cols){
#  data[,col] = as.integer(impute(data[,col],mean))
#}
