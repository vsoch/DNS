# This script will look at the different distributions of behavioral questions.  We simply want to see
# if we have normal, skewed, etc.

setwd("/home/vanessa/Documents/Work/DNS")
load("DNSBattery.Rda")
outpdf = "BehavioralDistributions.pdf"

# Get rid of the helper text
helpertext = c("Please use the scale to indicate how true each statement is.","This questionnaire has been designed to investigate ideas about intelligence. There are no right or wrong answers. We are interested in your ideas. Using the scale below, please indicate the extent...","Please select the answer that you most identify with:","Please answer the following questions or indicate whether you agree with the statements below.","Please answer the questions below, rating yourself on each of the criteria shown using the scale. As you answer each question, select the response that best describes how you have felt and conducte...","For each item, indicate how much you agree or disagree with what  the item says.  Respond to each...","As normal, healthy individuals go through their typical day they experience pain from time to time....","Please answer the following questions thinking about the most recent time you have had bodily pain s..","Please use the scale to indicate how accurately each statement describes how you GENERALLY ARE relat...","Please continue to use the rating scale next to each phrase to describe how accurately each statemen...","Please indicate if you have experienced any of the following events over the PAST YEAR. Do NOT select an event if it happened longer than 1 year ago.","Please continue to use the scale to select how frequently specific events occur or have occurred in...","Please select the frequency that you engage in each activity:","For each of the remaining questions, check the one best response. Please answer all questions.During...","Please indicate how much you agree with each statement.","Indicate how often you behave in the stated manner.","Please indicate how often you engaged in these behaviors in the past year. Please be honest in answe..","Please use the scale to select what is GENERALLY TRUE for you.","Please answer the following questions.","Please use the scale to rate how true each statement is for you.")
labels = battery$columns$QUESTION
for (text in helpertext){
  labels = gsub(text,"",labels)
}

battery$labels = labels
save(battery,file="DNSBattery.Rda")

# Visualize distributions
# Note that the label have some illegal characters - will be substituted automatically
pdf(outpdf,onefile=TRUE)
for (col in 1:ncol(battery$battery)){
  tmp = battery$battery[,col]
  tmp = as.numeric(as.character(tmp[which(tmp!=".")]))
  if (length(tmp)>0){
    label = battery$labels[match(colnames(battery$battery)[col],battery$columns$LABEL)]
    description = as.character(battery$labels[match(colnames(battery$battery)[col],battery$columns$LABEL)])
    # Get a random color FUN!
    color =  sample(colours(), 1)
    hist(tmp,col=color,main=label,sub=description,xlab="")
  }
}
dev.off()


# Sparse hierarchical clustering
tmp = as.data.frame(bat,stringsAsFactors=FALSE)
heatmap(tmp)

for (col in 1:ncol(tmp)){
  if (class(tmp[,col]) != "numeric"){
    tmp[,col] = as.numeric(tmp[,col])
    cat(col,"|")    
  }
}

battery$matrix = tmp
battery$README = c("filtered: raw data with . replaced with NA\nbattery: raw data no changes\nmatrix:numeric matrix\ntextidx: text or check questions")

# UNSUPERVISED CLUSTERING ---------------------------------------------

# Normalize values
norm = scale(battery$matrix)
battery$norm = norm


# We need to calculate distance matrix based on questions that both participants have answered
# We don't want to see clustering based on when they took the battery!
dissimilarity = array(dim=c(nrow(battery$matrix),nrow(battery$matrix)))
colnames(dissimilarity) = rownames(battery$matrix)
rownames(dissimilarity) = rownames(battery$matrix)
for (i in 1:nrow(battery$norm)){
  cat("Processing",i,"of",nrow(battery$norm),"\n")
  for (j in 1:nrow(battery$norm)){
    rowi = which(!is.na(battery$norm[i,]))
    rowj = which(!is.na(battery$norm[j,]))
    idx = intersect(rowi,rowj)
    dissimilarity[i,j] = dist(rbind(as.numeric(battery$norm[i,idx]),as.numeric(battery$norm[j,idx])))
  }
}
battery$dist = list()
battery$dist$euclidean = dissimilarity
save(battery,file="DNSBattery.Rda")

# Compare distance matrices
image(battery$dist$euclidean,main="Distance Matrix Accounting for Different Timeperiods",xaxt="n",yaxt="n")
image(as.matrix(dist(battery$matrix)),main="Distance Matrix Not Accounting for Different Timeperiods",xaxt="n",yaxt="n")

library(sparcl)
tmp = scale(battery$matrix)
# Do tuning parameter selection for sparse hierarchical clustering
perm.out <- HierarchicalSparseCluster.permute(tmp, wbounds=c(1.5,2:6),nperms=5)
print(perm.out)
pdf("sparsehcPlot.pdf")
plot(perm.out)
dev.off()

# Perform sparse hierarchical clustering
sparsehc = HierarchicalSparseCluster(dists=perm.out$dists,wbound=perm.out$bestw, method="complete")

sparsehc = HierarchicalSparseCluster(dists=as.matrix(disty), method="complete")
model = list(data=battery$matrix,sparsehc=sparsehc,disty=disty)
save(model,file="sparseHCmodel.Rda")

plot(sparsehc)
print(sparsehc)
str(hc)

# Plot using knowledge of class labels in order to compare true class
#   labels to clustering obtained
par(mfrow=c(1,1))
ColorDendrogram(sparsehc$hc,y=y,main="My Simulated Data",branchlength=.007)
# Now, what if we want to see if out data contains a *secondary*
#   clustering after accounting for the first one obtained. We
#   look for a complementary sparse clustering:
sparsehc.comp <- HierarchicalSparseCluster(x,wbound=perm.out$bestw,
                                           method="complete",uorth=sparsehc$u)
# Redo the analysis, but this time use "absolute value" dissimilarity:
perm.out <- HierarchicalSparseCluster.permute(x, wbounds=c(1.5,2:6),
                                              nperms=5, dissimilarity="absolute.value")
print(perm.out)
plot(perm.out)
# Perform sparse hierarchical clustering
sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists, wbound=perm.out$bestw, method="complete", dissimilarity="absolute.value")
par(mfrow=c(1,2))
plot(sparsehc)

# TRY THIS PACKAGE FOR CLUSTERING
library("RGLUEANN")
http://rogiersbart.blogspot.com/2014/10/rglueann-package-available-on-github.html
demo("RGLUEANN_training_and_prediction")
demo("RGLUEANN_cross-validation")
