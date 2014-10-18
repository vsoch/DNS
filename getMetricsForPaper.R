setwd("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults")
setwd("C:\\Users\\Vanessa\\Documents\\Work\\GENE_EXPRESSION\\neurosynth\\sigresult")

# Numbe of unique samples in ABA
samp = read.csv("C:\\Users\\Vanessa\\Documents\\Work\\ALLEN\\allen-brain-Samples.csv",head=TRUE)
length(unique(samp$structure_name))
# 414

# Number of significant datasets per disorder
load("disorderGeneListsN8pt1.Rda")
disorder = unique(filter$DISORDER)
for (dd in 1:length(disorder)){
  d=disorder[dd]
  subset = filter[filter$DISORDER == d,]
  num = unique(subset$FOLDER)
  cat(as.character(d),"\n")
  cat(as.character(num))
  cat("\n\n")
}

# ALZ 
# GSE12685_ADvsHC_brainTerms_genes.Gsea.1405202376693 
# GSE1297_incipADvsHC_brainTerms.Gsea.1402269722957 
# GSE16759GPL570_parietalADvsHC_brainTerms.Gsea.1402336352587 
# GSE29378_hippoADvsHC_brainTerms.Gsea.1402184255413

# ASD 
# GSE38322_ASDvsHC_cer_brainTerms.Gsea.1402181385337

# LUPUS 
# GSE22098LUPUSvsHC_brainTerms.Gsea.1404769742065 
# GSE30153LUPUSvsHC_brainTerms.Gsea.1404770822137

# MS 
# GSE17048MSvsHC_brainTerms.Gsea.1404845552954 
# GSE41849MSvsHC_brainTerms.Gsea.1404845551566

# PARK 
# GSE28894PARKCERvsHC_brainTerms.Gsea.1405025664761

# PTSD 
# GSE860_PTSDvsHC_brainTerms.Gsea.1402256689658

# RETT 
# GSE34099RETTvsHC_brainTerms.Gsea.1405025789769

# SZO 
# GSE12679_SZOvsHC_brainTerms.Gsea.1402263848076

# Number of genes core for each disorder (Table 1)
load("disorderGenespt1.Rda")
coreGeneCounts = c()
for (g in unique(genes$DISORDER)){
  coreGeneCounts[g] = length(unique(genes$PROBE[which(genes$DISORDER == g)]))
  cat(g,":",coreGeneCounts[g],"\n")
}
save(coreGeneCounts,file="coreGeneCountspt1.Rda")

# How many gene probes in entire Allen Brain Atlas data?
load("/home/vanessa/Documents/Work/ALLEN/GeneExpressionMeans3702.Rda")
# how many genes?
dim(expressionMeans)
load("/home/vanessa/Documents/Work/ALLEN/probesExpressionAll.Rda")
# 29131 gene symbols
# How many probes?
dim(probes)
# 58871

# The 525 brainterm maps ranged in size from N to N sample points
setwd("/home/vanessa/Documents/Work/GEgNE_EXPRESSION/neurosynth/probeSets/9mmsq/data")
files = list.files(pattern="_up.Rda")
maxn = 0
minn = 999
sizes = c()
for (f in files){
  load(f)
  minn = min(nrow(mat),minn)
  maxn = max(nrow(mat),maxn)  
  sizes = c(sizes,nrow(mat))
}
# 27 to 952
save(sizes,file="/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/brainTermSetSizes.Rda")
# Of a potential 1050 brainterm gene sets there were this many significant ranging in size from N to N

# We systematically compared each of the this many brainterm gene sets to the pattern of gene expression in all 56 neurological disorders
# It is notable that a majority of our significant results (HOW MANY?) 

genesets = "C:\\Users\\Vanessa\\Documents\\Dropbox\\Code\\R\\GSEA-P-R\\GeneSetDatabases\\brainTerms.gmt"

x <- scan(genesets, what="", sep="\n")
# Separate elements by one or more whitepace
y <- strsplit(x, "\t")
# Extract the first vector element and set it as the list element name
names(y) <- sapply(y, `[[`, 1)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
y <- lapply(y, `[`, -1)
#y <- lapply(y, function(x) x[-1]) # same as above
note = "Allen Brain Atlas expression up and down regulated subsets for neurosynth fdr .05 corrected brain map sorted by pval more sig first"
n2g = list(genesets=y,note=note)
save(n2g,file="n2g_genesets.Rda")

# What are the different gene set sizes?
counts = c()
for (t in 1:length(n2g$genesets)){
  counts = c(counts,length(n2g$genesets[[t]]))  
}

# Let's get the significant GENE set data input for each:
list.files()
genesets = unique(genes$FOLDER)
cat(as.character(genesets),file="DataSetsSignicantpt1.txt",sep="\n")

# Load genes for all 8 disorders
load("disorderGenesN8pt1.Rda")
disorders = as.character(unique(genes$DISORDER))
for (d in disorders){
  subset = genes[which(genes$DISORDER==d),]
  g = names(sort(table(as.character(subset$PROBE)),decreasing=TRUE)) # sorted counts
  c = as.numeric(sort(table(as.character(subset$PROBE)),decreasing=TRUE))
  cat(g,sep="\n")
  cat(c,sep="\n")
}

write.table(genes,row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t",file="disorderGenespt1.tab")
getwd()

# Number of significant terms for each disorder - plus som
# ALZ UP: 304 DOWN: 0 
# ASD UP: 2 DOWN: 0 
# LUPUS UP: 25 DOWN: 0 
# MS UP: 3 DOWN: 196 
# PARK UP: 1 DOWN: 0 
# PTSD UP: 1 DOWN: 0 
# RETT UP: 3 DOWN: 0 
# SZO UP: 1 DOWN: 0
setwd("C:\\Users\\Vanessa\\Documents\\Dropbox\\Research\\Share\\Diabolical\\gsea\\brainterms\\MonsterGSEA")
load("disorderGeneListsN8pt1.Rda")
disorder = as.character(unique(filter$DISORDER))
for (d in disorder){
  subset = filter[filter$DISORDER == d,]
  up = gsub("_UP","",unique(as.character(subset$TERM[grep("_UP",subset$TERM)])))
  down = gsub("_DOWN","",unique(as.character(subset$TERM[grep("_DOWN",subset$TERM)])))
  cat(as.character(d),as.character("UP:"),length(up),as.character("DOWN:"),length(down),"\n")
}

# Ready SOM so we can plot terms on it
# This script will create an (R) som image by reading in
# match score files in /scratch/PI/dpwall/DATA/IMAGING_GENOMICS/BRAINGRID/score/
library(RColorBrewer)
setwd("C:\\Users\\Vanessa\\Documents\\Dropbox\\Code\\Google\\SGE\\braingrid\\maps")
load("brainMap.Rda")

# We need to define a color scale that indicates the strength of the match score
colorscale = brewer.pal(9,"YlOrRd")
colorscale = colorRampPalette(brewer.pal(8,"YlOrRd"))(100)

# Prepare index of where each brainterm is on the map
lookup = list()
for (l in 1:length(brainMap$labels)){
  label = strsplit(brainMap$labels[l],"\n")[[1]]
  for (i in label){
    lookup[i] = l
  }  
}

# Create three color gradient, because we have 3 levels of disorder significance
colorscale = c("#F5DA81","#FAAC58","#FF8000")

# For each disorder, plot where the term maps are
for (d in disorder){
  subset = filter[filter$DISORDER == d,]
  up = gsub("_UP","",(as.character(subset$TERM[grep("_UP",subset$TERM)])))
  down = gsub("_DOWN","",(as.character(subset$TERM[grep("_DOWN",subset$TERM)])))
  all = table(c(up,down))

  # Create list of relevant labels
  colors = array("#FFFFFF",dim=nrow(brainMap$som$grid$pts))
  for (aa in 1:length(all)){
    a = all[aa]
    colors[as.numeric(lookup[tolower(names(a))])] = colorscale[as.numeric(a)]
  }
  png(file=paste(d,"_map.png",sep=""),width=14,height=14,units="in",res=300)
  plot(brainMap$som$grid$pts,main=paste("Significant BrainTerm Maps for",d),col=colors,xlab="Nodes",ylab="Nodes",pch=15,cex=8)
  text(brainMap$som$grid$pts,brainMap$labels,cex=.6)
  dev.off()
}


library("Rniftilib")
library("RColorBrewer")

# Define distance functions
euc.dist = function(x1,x2) {sqrt(sum((x1 - x2) ^ 2))}
cos.dist = function(x1,x2) {crossprod(x1, x2)/sqrt(crossprod(x1) * crossprod(x2))}
allscores=c()
# Now let's read in all the maps for a disorder, take a mean, and plot in space
load("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/disorderGeneListsN8pt1.Rda")
disorder = as.character(unique(filter$DISORDER))
load("/home/vanessa/Documents/Dropbox/Code/Google/SGE/braingrid/R/brainMap.Rda")
mapdir = "/home/vanessa/Documents/Work/BRAINSPAN/dim_reduction/nsynth525"
setwd(mapdir)
maps = list.files(mapdir,pattern="8mm*")
# For each disorder, plot where the term maps are
for (d in disorder){
  subset = filter[filter$DISORDER == d,]
  up = gsub("_UP","",(as.character(subset$TERM[grep("_UP",subset$TERM)])))
  down = gsub("_DOWN","",(as.character(subset$TERM[grep("_DOWN",subset$TERM)])))
  all = tolower(unique(c(up,down)))
  alltable = table(c(up,down))
  
  # Read in each image to matrix
  data = matrix(nrow=length(all),ncol=16128)
  for (a in 1:length(all)){
    nii = nifti.image.read(paste("8mm",all[a],".nii",sep=""),read_data=1)
    niivector = as.vector(nii[,,,1])
    data[a,] = niivector
  }
  
  # Take mean image
  newnii = colMeans(data)
  
  # Convert to z score, and try taking only op 5% of values (Z=1.96)
  #Z = (newnii - mean(newnii[newnii!=0])) / sd(newnii[newnii!=0])
  #Z[abs(Z)<1.96] = 0
  template = nii
  template[,,] = niivector
  nifti.set.filenames(template,paste("/home/vanessa/Desktop/",d,"_mean.nii",sep=""))
  nifti.image.write(template)  

  # Save data to file, in case we want to use it later
  #brainTerms = list(data=data,labels=all,counts=alltable,mean=newnii,Z=Z)
  #save(brainTerms,file=paste(d,"_allBrainTerms.Rda",sep=""))
  
  # Now to create som - match mean image to each node in SOM
  matchscores = array(0,ncol(brainMap$som$data))
  
  for (i in 1:nrow(brainMap$som$codes)){
    node = brainMap$som$codes[i,]
    matchscores[i] = euc.dist(node,newnii)
  }
  # Create list of relevant labels
  rbPal <- colorRampPalette(brewer.pal(8,"YlOrRd"))
  colors = rbPal(10)[as.numeric(cut(matchscores,breaks = 10))]
  
  png(file=paste("/home/vanessa/Desktop/",d,"_meanTermsZ05.png",sep=""),width=14,height=14,units="in",res=300)
  plot(brainMap$som$grid$pts,main=paste("Mean of Significant BrainTerm Maps for",d),col=colors,xlab="Nodes",ylab="Nodes",pch=15,cex=8)
  text(brainMap$som$grid$pts,brainMap$labels,cex=.6)
  dev.off()
}


# Try clustering mean images for each disorder
setwd("C:\\Users\\Vanessa\\Documents\\Work\\GENE_EXPRESSION\\neurosynth\\sigresult\\visualization\\meanTermImages")
files = list.files(pattern="*.Rda")
meanys = array(dim=c(8,16128))
for (f in 1:length(files)){
  ff = files[f]
  load(ff)
  meanys[f,] = brainTerms$mean
}

disty = dist(meanys)
hc = hclust(disty)
plot(hc)

# SPARSE HIERARCHICAL CLUSTERING TO DO FEATURE SELECTION
setwd("/scratch/users/vsochat/DATA/GENE_EXPRESSION/neurosynth/sigresult/disorderBrains")
load("/scratch/users/vsochat/DATA/BRAINMAP/dimensionality_reduction/som_pAgF/brainMap.Rda")
load("SZO_allBrainTerms.Rda")
cos.dist = function(x1,x2) {crossprod(x1, x2)/sqrt(crossprod(x1) * crossprod(x2))}

library("sparcl")
library("RColorBrewer")
library("Rniftilib")
#idx = which(colSums(brainTerms$data)!=0)
#subset = brainTerms$data[,idx]
res = HierarchicalSparseCluster(x=brainTerms$data,niter=5)
# Make into z scores, get top .05 of weights
test = c(res$ws,-res$ws)
#test = brainTerms$data
Z = (test - mean(test[test!=0])) / sd(test[test!=0])
idx = which(Z>2.4)
# ALZ: [1] 0.04347603 to [1] 0.1821919
# MS: [1] 0.0408482 to [1] 0.1944342
# LUPUS: [1] [1] 0.06724099 to [1] 0.2124244
# ASD: ] 0.1801578 to 0.3575494
# RETT:  0.1207415 to [1] 0.4273447
# PARK AND SZO are single vectors
save(res,file="PARK_sparcl.Rda")
# Create a new image of same size, take mean over important voxels
niivector = array(0,dim=ncol(brainTerms$data))
# Read in each image to matrix
library("Rniftilib")
nii = nifti.image.read("/scratch/users/vsochat/DATA/BRAINMAP/nsynth525pFgA/8mm/r1back.nii",read_data=1)
template = nii
#niivector[idx] = brainTerms$data[idx]
niivector[idx] = colMeans(brainTerms$data[,idx])
template[,,] = niivector
nifti.set.filenames(template,"SZO_sparseHCbrain2pt4.nii")
nifti.image.write(template)  
# Create som
matchscores = array(0,nrow(brainMap$som$codes))
for (i in 1:nrow(brainMap$som$codes)){
  node = brainMap$som$codes[i,]
  matchscores[i] = cos.dist(node,niivector)
}
matchscores = c(0,matchscores,1)
# Create list of relevant labels
rbPal <- colorRampPalette(brewer.pal(8,"YlOrRd"))
colors = rbPal(10)[as.numeric(cut(matchscores,breaks = 10))]
colors= colors[-c(1,508)]

png(file=paste("SZO_sparseTerms2pt4.png",sep=""),width=14,height=14,units="in",res=300)
plot(brainMap$som$grid$pts,main=paste("Summary of Significant BrainTerm Maps for SZO"),col=colors,xlab="Nodes",ylab="Nodes",pch=15,cex=8)
text(brainMap$som$grid$pts,brainMap$labels,cex=.6)
dev.off()
