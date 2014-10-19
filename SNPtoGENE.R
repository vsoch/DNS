# This script will read in a matrix of ndar (or DNS) SNP data,
# map those SNPs to genes, and cross list the genes
# with the neuro2gene brainterm genes.  We
# will save a matrix for each subset of genes
# to do an analysis - eventually looking for
# correlation between these genes and behavioral
# metrics for this group

setwd("/scratch/users/vsochat/DATA/GENE_EXPRESSION/ndar")
setwd("/home/vanessa/Documents/Work/DNS")
files = list.files(pattern="*.csv")
snps = c()
for (f in files){
  raw = read.table(f,sep=",",head=FALSE,skip=1)
  raw = raw[,c(1,2)]
  raw = raw[-which(raw$V2==""),]
  snps = rbind(snps,raw)
}

# Save list of SNPs
colnames(snps) = c("SNP","GENE")
save(snps,file="rawSNPList.Rda")

# Now create a dictionary with gene lookup
snp = unique(snps$SNP)
snp.lookup = list()
for (s in 1:length(snp)){
  snp.lookup[[snp[s]]] =
  
}

# Get all the snp names to look up
lookup = as.character(raw$ID)

# Now read in my list of genes
load("/home/vanessa/Documents/Work/GENE_EXPRESSION/neurosynth/results/sigresults/disorderGenesN8pt1.Rda")
dnsgenes = unique(as.character(snps$GENE))
mygenes = unique(as.character(genes$PROBE))
mygenes[which(mygenes %in% dnsgenes)]

overlap = mygenes[which(mygenes %in% dnsgenes)]

# filter down DNS genes to my genes
dns = snps[which(as.character(snps$GENE) %in% overlap),]
save(dns,file="dnsGeneOverlap.Rda")

# Woohoo!

# Filter down to autism genes (or all)
asd_genes = as.character(uniquegenes$ASD)
uniquegenes = uniquegenes[c(3,4,10,12,13,15,17,18)]
genes = c()
for (l in 1:length(uniquegenes)){
  genes = c(genes,as.character(uniquegenes[[l]]))
}
genes = unique(genes)

# Find overlap
idx = which(mapping$gene %in% genes)

# For ALL genes, we have 52.  Let's create a matrix of just
# Those genes
keepers = mapping[idx,]
filter = raw[which(raw$ID %in% keepers$snp),]
snp = list(mapping=keepers,data=filter)
save(snp,file="n2g_SNP_matrix.Rda")

# STOPPED HERE - need to parse behavioral data!
