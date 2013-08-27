#install.packages("mvnormtest")
#install.packages("npmv")
list.files()
library(mvnormtest)
library(imputation)
library(npmv)


arc<-read.csv("ARC TLDA for stats 7-20-13.csv") # load the data
poa<-read.csv("mPOA TLDA for stats 7-20-13.csv")

head(arc)

attach(arc)
names(arc)
which(colnames(arc)=="E2.RIA") #start of hormone data
which(colnames(arc)=="ACTH") #end
hormone<-arc[,c(2:8,66:73)]

hormone$Group<-unlist(lapply(list(hormone$Group),FUN=gsub,pattern=" ",replacement="_"))

shapiro.test(log10(E2.RIA))
shapiro.test(log10(LH))
shapiro.test(log10(FSH))
shapiro.test(log10(GH))
shapiro.test(sqrt(TSH))
shapiro.test(log10(PRL))
shapiro.test(log10(BDNF))
shapiro.test(log10(ACTH))
detach(arc)

cols<-which(is.na(hormone),arr.ind=TRUE)

imputed<-kNNImpute(log10(hormone[,8:15]),3)$x
mshapiro.test(t(as.matrix((imputed)))) # not normal or multivariate normal

imputed$Group<-hormone$Group

#nonparametric anova
nonpartest(E2.RIA|LH|FSH|GH|TSH|PRL|BDNF|ACTH~Group,data=imputed)
hack1<-hackNonParTest(E2.RIA|LH|FSH|GH|TSH|PRL|BDNF|ACTH~Group,data=imputed,permreps
                      =50000)
wtf<-ssnonpartest(E2.RIA|LH|FSH|GH|TSH|PRL|BDNF|ACTH~Group,data=imputed)


man1<-manova(cbind(E2.RIA,LH,FSH,GH,TSH,PRL,BDNF,ACTH)~Group)#not normal
summary(man1,test="Pi")
summary.aov(man1)
ACTHanova<-aov(ACTH~Group)
Pairs <- glht(ACTHanova, linfct = mcp(Group = "Tukey"))
cld(Pairs)


library(vegan)
install.packages('imputation')

pman1<-adonis(imputed$x~Group)

install.packages("MCMCglmm")

library(MASS)
noSwitch<-arc[1:60,]

fitJack<-lda(Group~E2.RIA+LH+FSH+GH+TSH+PRL+BDNF+ACTH,na.action="na.omit",CV=TRUE,data=data.frame(imputed$x))
ct <- table(Group, fitJack$class)

fitNo<-lda(Group~E2.RIA+LH+FSH+GH+TSH+PRL+BDNF+ACTH,na.action="na.omit",data=(noSwitch))

fit<-lda(Group~E2.RIA+LH+FSH+GH+TSH+PRL+BDNF+ACTH,na.action="na.omit",data=data.frame(scale(cbind(E2.RIA+LH+FSH+GH+TSH+PRL+BDNF+ACTH))))
plot(fit,dimen=2)