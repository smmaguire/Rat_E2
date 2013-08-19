#GeneExpression

#############TO_DO########
# make hclust heat maps to show the fact that the genes cluster by brain area rather than be treatment so much
#
#
##########################
#Covariance networks between the two brain areas
arc<-read.csv("ARC TLDA for stats 7-20-13.csv") # load the data
poa<-read.csv("mPOA TLDA for stats 7-20-13.csv")

names(arc)[11]
names(arc)[57] #bounds are 11 to 57

arc.mat<-arc[,11:57] # column 11 is missing data.
which(grepl("Mtnr1a",names(arc.mat)))
arc.mat[,30]
which(arc.mat=="#VALUE!",arr.ind=TRUE)
arc.mat[74,30]<-NA
arc.mat[,30]<-as.numeric(as.character(arc.mat[,30]))
arc.mat.small<-(arc.mat[,c(-11,-31,-33)]) # this is to make it match witht he POA data
arc.mat<-(arc.mat[,-11])
arc.mat<-as.matrix(arc.mat)

poa.mat<-(poa[,11:57])
str(poa.mat) #Dbh, Mtnr1b, Npvf
which(grepl("Dbh",names(poa.mat))) #11
which(grepl("Mtnr1b",names(poa.mat))) #31
which(grepl("Npvf",names(poa.mat))) #33

poa.mat<-as.matrix(poa.mat[,c(-11,-31,-33)])


poa.cov<-cov(poa.mat,use="pairwise.complete.obs")
arc.cov<-cov(arc.mat.small,use="pairwise.complete.obs")

cor(as.numeric(diag.remove(poa.cov)),as.numeric(diag.remove(arc.cov)),use="pairwise.complete.obs")

#should probably consider E2 vs non-E2 treatment + brain area
arc.mat.e2.1<-arc.mat.small[arc$E2.effect==1,]
arc.mat.e2.2<-arc.mat.small[arc$E2.effect==2,]

poa.mat.e2.1<-poa.mat[poa$E2.effect==1,]
poa.mat.e2.2<-poa.mat[poa$E2.effect==2,]

arc.mat.e2.1.cov<-cov(arc.mat.e2.1,use="pairwise.complete.obs")
poa.mat.e2.1.cov<-cov(poa.mat.e2.1,use="pairwise.complete.obs")

arc.mat.e2.2.cov<-cov(arc.mat.e2.2,use="pairwise.complete.obs")
poa.mat.e2.2.cov<-cov(poa.mat.e2.2,use="pairwise.complete.obs")

e2.1.across_areas<-cor(as.numeric(diag.remove(arc.mat.e2.1.cov)),as.numeric(diag.remove(poa.mat.e2.1.cov)),use="pairwise.complete.obs")
e2.2.across_areas<-cor(as.numeric(diag.remove(arc.mat.e2.2.cov)),as.numeric(diag.remove(poa.mat.e2.2.cov)),use="pairwise.complete.obs")

arc.across_treatment<-cor(as.numeric(diag.remove(arc.mat.e2.1.cov)),as.numeric(diag.remove(arc.mat.e2.2.cov)),use="pairwise.complete.obs")
poa.across_treatment<-cor(as.numeric(diag.remove(poa.mat.e2.1.cov)),as.numeric(diag.remove(poa.mat.e2.2.cov)),use="pairwise.complete.obs")

###################################################WGCNA##################################################################

install.packages('WGCNA')
library(WGCNA)
allowWGCNAThreads() # faster w/ multiple threads

#1 fix names.
gene<-vector()
for(probe in colnames(arc.mat)){
  gsplt<-unlist(strsplit(probe,split=".",fixed=TRUE))
  gene<-c(gene,gsplt[1])
}

colnames(arc.mat)<-gene
rownames(arc.mat)<-as.character(arc$Sample.Name)
gsg <- goodSamplesGenes(arc.mat, verbose = 3);

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(arc.mat)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(arc.mat)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  arc.mat.gsg <- arc.mat[gsg$goodSamples, gsg$goodGenes]
}
#Removing samples: ABV26, ABV75, ABV49, ABV44, ABV45, ABV46, ABV55

sampleTree = flashClust(dist(arc.mat.gsg), method = "average");
sizeGrWindow(12,9)
#Removing ABV86
arc.mat.gsg<-arc.mat.gsg[rownames(arc.mat.gsg)!="ABV86",]
nGenes<-ncol(arc.mat.gsg)
nSamples<-nrow(arc.mat.gsg)

sampNames<-rownames(arc.mat.gsg)
traitRows<- match(sampNames,as.character(arc$Sample.Name))
nas<-which(is.na(arc$E2.effect))
arc$E2.effect[nas]<-3
arc$age.effect
datTraits<- data.frame(Treatment=as.numeric(arc$Group[traitRows]),E2=arc$E2.effect[traitRows],age=arc$age.effect[traitRows])

sampleTree2 = flashClust(dist(arc.mat.gsg), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE);

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")