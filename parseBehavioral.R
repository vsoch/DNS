# This script will parse behavioral data from the Duke NeuroGenetics Study (DNS)
# We will read in the data from the xlsx provided by Annchen and Spenser, organize into
# RDA objects, and then visualize the distributions.

library("xlsx")
setwd("/home/vanessa/Documents/Work/DNS")

# DATA DICTIONARY ----------------------------------------------------------------------

# CLINICAL DATA ------------------------------------------------------------------------
clinical = "Vanessa_Day1.xlsx"
neuropsych = read.xlsx(clinical, 1)
diagnoses = read.xlsx(clinical, 2)
meds = read.xlsx(clinical, 3)
sdmt = read.xlsx(clinical, 4)
clinical = list(neuropsych=neuropsych,diagnoses=diagnoses,meds=meds,sdmt=sdmt)

# Find empty columns and get rid of
for (m in 1:length(clinical)){
  datatype = names(clinical)[m]
  tmp = clinical[[datatype]]
  todelete = c()
  for (c in ncol(tmp)){
     if(all(is.na(neuropsych[,c]))){
       todelete = c(todelete, c)
       cat(datatype,": Adding column",names(neuropsych[,c]),"to delete list, column #",c,"\n")
     }
  }
  # Do the deletion
  if (length(todelete)>0){
    for (col in todelete){
      clinical[[datatype]] = clinical[[datatype]][,-col]
    }
  }
}

# Actually, it's probably better to do this manually
clinical$neuropsych = clinical$neuropsych[,-seq(48,57)]
clinical$diagnoses = clinical$diagnoses[,-41]

# Set row name to be dns id
rownames(clinical$neuropsych) = clinical$neuropsych$dns_id
rownames(clinical$meds) = clinical$meds$dns_id
rownames(clinical$diagnoses) = clinical$diagnoses$dns_id
rownames(clinical$sdmt) = clinical$sdmt$dns_id

clinical$neuropsych = clinical$neuropsych[,-1]
clinical$meds = clinical$meds[,-1]
clinical$diagnoses = clinical$diagnoses[,-1]
clinical$sdmt = clinical$sdmt[,-1]

# save to file
save(clinical,file="DNSClinical.Rda")

# BEHAVIORAL DATA ----------------------------------------------------------------------

# The data is too big to read from xlsx, save to tab separated file
battery = read.csv("battery.tab",sep="\t")
rownames(battery) = battery$ID
battery = battery[,-1]
labels = read.csv("batteryGroups.tab",sep="\t")
#idx = which(labels$LABEL %in% colnames(battery))
#labels = labels[idx,]
battery = list(battery=battery,columns=labels)
save(battery,file="DNSBattery.Rda")

# Add filtered data - has "." replaced with NA
filtered = read.csv("batteryfilter.tab",sep="\t")
rownames(filtered) = filtered$ID
filtered = filtered[,-1]
battery$filtered = filtered
save(battery,file="DNSBattery.Rda")
