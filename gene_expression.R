#GeneExpression
#wtf
#############TO_DO########
# make hclust heat maps to show the fact that the genes cluster by brain area rather than be treatment so much
# bagged regression to try and predict the different effects/ or some sort of machine learning
# 
##########################
#Covariance networks between the two brain areas
library(sna)
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
source("http://bioconductor.org/biocLite.R")
biocLite("impute")
install.packages("WGCNA") 
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

powers<- c(1:10, seq(from=12,to=20,by=2))
sft<-pickSoftThreshold(arc.mat.gsg,powerVector=powers,verbose=5)
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1=.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="soft threshold (power)",ylab="scale free R2",type="n",main="scale independence")
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=.9,col="red")

plot(sft$fitIndices[,1],sft$fitIndices[,5],type="n",xlab="power",ylab="mean connectivity")
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,cex=cex1,col="red")

#I'll go with 4 for now but may test from 4-8

softPower<-5
adjacency=adjacency(arc.mat.gsg,power=softPower)
TOM=TOMsimilarity(adjacency)
dissTOM<- 1-TOM

geneTree=flashClust(as.dist(dissTOM),method="average")
sizeGrWindow(12,9)
plot(geneTree,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity",labels=FALSE,hang=.04)

minModuleSize=0
dynamicMods<- cutreeDynamic(dendro=geneTree,distM=dissTOM,deepSplit=2,pamRespectsDendro=FALSE,minClusterSize=minModuleSize)
table(dynamicMods)

dynamicColors=labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree,dynamicColors,"Dynamic Tree Cut",dendroLabels=FALSE,hang=.03,addGuide=TRUE,guideHang=.05,main="Gene dendrogram and module colors")

MEList<-moduleEigengenes(arc.mat.gsg,colors=dynamicColors)
Modules<-data.frame(gene=colnames(arc.mat.gsg),color=MEList$validColors)
MEs<-MEList$eigengenes

MEs$treatment<-arc$Group[match(rownames(arc.mat.gsg),arc$Sample.Name)]
blue<-aov(MEblue~treatment,data=MEs)
summary(blue)
plot(MEblue~treatment,data=MEs)

brown<-aov(MEbrown~treatment,data=MEs)
summary(brown)
plot(MEbrown~treatment,data=MEs)

green<-aov(MEgreen~treatment,data=MEs)
summary(green)
plot(MEgreen~treatment,data=MEs)

names(MEs)
grey<-aov(MEgrey~treatment,data=MEs)
summary(grey)

red<-aov(MEred~treatment,data=MEs)
summary(red)
plot(MEred~treatment,data=MEs)

turquoise<-aov(MEturquoise~treatment,data=MEs)
summary(turquoise)
plot(MEturquoise~treatment,data=MEs)

yellow<-aov(MEyellow~treatment,data=MEs)
summary(yellow)
plot(MEyellow~treatment,data=MEs)


##################################POA###############################
#1 fix names.
#rm(list=ls())
gene<-vector()
for(probe in colnames(poa.mat)){
  gsplt<-unlist(strsplit(probe,split=".",fixed=TRUE))
  gene<-c(gene,gsplt[1])
}

colnames(poa.mat)<-gene
rownames(poa.mat)<-as.character(poa$Sample.Name)
gsg <- goodSamplesGenes(poa.mat, verbose = 3);

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(poa.mat)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(poa.mat)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  poa.mat.gsg <- poa.mat[gsg$goodSamples, gsg$goodGenes]
}
#Removing samples: ABV26, ABV85, ABV18, ABV49, ABV43, ABV56

sampleTree = flashClust(dist(poa.mat.gsg), method = "average");
sizeGrWindow(12,9)
plot(sampleTree)

nGenes<-ncol(poa.mat.gsg)
nSamples<-nrow(poa.mat.gsg)

sampNames<-rownames(poa.mat.gsg)
traitRows<- match(sampNames,as.character(poa$Sample.Name))
nas<-which(is.na(poa$E2.effect))
poa$E2.effect[nas]<-3
poa$age.effect
datTraits<- data.frame(Treatment=as.numeric(poa$Group[traitRows]),E2=poa$E2.effect[traitRows],age=poa$age.effect[traitRows])

sampleTree2 = flashClust(dist(poa.mat.gsg), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE);

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

powers<- c(1:10, seq(from=12,to=20,by=2))
sft<-pickSoftThreshold(poa.mat.gsg,powerVector=powers,verbose=5)
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1=.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="soft threshold (power)",ylab="scale free R2",type="n",main="scale independence")
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=.9,col="red")

plot(sft$fitIndices[,1],sft$fitIndices[,5],type="n",xlab="power",ylab="mean connectivity")
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,cex=cex1,col="red")

#I'll go with 4 for now but may test from 4-8

softPower<-5
adjacency=adjacency(poa.mat.gsg,power=softPower)
TOM=TOMsimilarity(adjacency)
dissTOM<- 1-TOM

geneTree=flashClust(as.dist(dissTOM),method="average")
sizeGrWindow(12,9)
plot(geneTree,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity",labels=FALSE,hang=.04)

minModuleSize=0
dynamicMods<- cutreeDynamic(dendro=geneTree,distM=dissTOM,deepSplit=2,pamRespectsDendro=FALSE,minClusterSize=minModuleSize)
table(dynamicMods)

dynamicColors=labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree,dynamicColors,"Dynamic Tree Cut",dendroLabels=FALSE,hang=.03,addGuide=TRUE,guideHang=.05,main="Gene dendrogram and module colors")

MEList<-moduleEigengenes(poa.mat.gsg,colors=dynamicColors)
Modules<-data.frame(gene=colnames(poa.mat.gsg),color=MEList$validColors)
MEs<-MEList$eigengenes

MEs$treatment<-poa$Group[match(rownames(poa.mat.gsg),poa$Sample.Name)]

names(MEs)
blue<-aov(MEblue~treatment,data=MEs)
summary(blue)
plot(MEblue~treatment,data=MEs)

brown<-aov(MEbrown~treatment,data=MEs)
summary(brown)

green<-aov(MEgreen~treatment,data=MEs)
summary(green)
plot(MEgreen~treatment,data=MEs)

grey<-aov(MEgrey~treatment,data=MEs)
summary(grey)
plot(MEgrey~treatment,data=MEs)

turquoise<-aov(MEturquoise~treatment,data=MEs)
summary(turquoise)

yellow<-aov(MEyellow~treatment,data=MEs)
summary(yellow)
plot(MEyellow~treatment,data=MEs)

black<-aov(MEblack~treatment,data=MEs)
summary(black)
plot(MEblack~treatment,data=MEs)

  magenta<-aov(MEmagenta~treatment,data=MEs)
summary(magenta)
plot(MEmagenta~treatment,data=MEs)

pink<-aov(MEpink~treatment,data=MEs)
summary(pink)

purple<-aov(MEpurple~treatment,data=MEs)
summary(purple)
plot(MEblack~treatment,data=MEs)

red<-aov(MEred~treatment,data=MEs)
summary(red)