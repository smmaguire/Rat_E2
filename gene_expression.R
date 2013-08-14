#GeneExpression

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
arc.mat<-(arc.mat[,-11])
arc.mat<-as.matrix(arc.mat)

poa.mat<-(poa[,11:57])
str(poa.mat) #Dbh, Mtnr1b, Npvf
which(grepl("Dbh",names(poa.mat))) #11
which(grepl("Mtnr1b",names(poa.mat))) #31
which(grepl("Npvf",names(poa.mat))) #33

poa.mat<-as.matrix(poa.mat[,c(-11,-31,-33)])

cov(poa.mat,use="pairwise.complete.obs")
