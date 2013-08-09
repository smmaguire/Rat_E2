setwd("/home/sean/Downloads/")
list.files()

arc<-read.csv("ARC TLDA for stats 7-20-13.csv") # load the data
poa<-read.csv("mPOA TLDA for stats 7-20-13.csv")

head(arc)

attach(arc)
install.packages("mvnormtest")
shapiro.test(log10(E2.RIA))
shapiro.test(log10(LH))
shapiro.test(log10(FSH))
shapiro.test(log10(GH))
shapiro.test(sqrt(TSH))
shapiro.test(log10(PRL))
shapiro.test(log10(BDNF))
shapiro.test(log10(ACTH))

cols<-which(is.na(hormone),arr.ind=TRUE)
hormone2<-hormone[,-as.vector(cols[,2])]
hormone2[c(1,2,3,4,6,7,8),]<-log10(hormone2[c(1,2,3,4,6,7,8),])
hormone2[5,]<-sqrt(hormone2[5,])
mshapiro.test(log10(hormone))

man1<-manova(cbind(E2.RIA,LH,FSH,GH,TSH,PRL,BDNF,ACTH)~Group)
summary(man1,test="Pi")
summary.aov(man1)
ACTHanova<-aov(ACTH~Group)
Pairs <- glht(ACTHanova, linfct = mcp(Group = "Tukey"))
cld(Pairs)


library(vegan)
install.packages('imputation')
library(imputation)
imputed<-kNNImpute(cbind(E2.RIA,LH,FSH,GH,TSH,PRL,BDNF,ACTH),3)
pman1<-adonis(imputed$x~Group)

install.packages("MCMCglmm")

library(MASS)
noSwitch<-arc[1:60,]

fitJack<-lda(Group~E2.RIA+LH+FSH+GH+TSH+PRL+BDNF+ACTH,na.action="na.omit",CV=TRUE,data=data.frame(imputed$x))
ct <- table(Group, fitJack$class)

fitNo<-lda(Group~E2.RIA+LH+FSH+GH+TSH+PRL+BDNF+ACTH,na.action="na.omit",data=(noSwitch))

fit<-lda(Group~E2.RIA+LH+FSH+GH+TSH+PRL+BDNF+ACTH,na.action="na.omit",data=data.frame(scale(cbind(E2.RIA+LH+FSH+GH+TSH+PRL+BDNF+ACTH))))
plot(fit,dimen=2)

anvM1<-aov(E2.RIA~Group)
Pairs <- glht(anvM1, linfct = mcp(Group = "Tukey"))
cld(Pairs)

shapiro.test(sqrt(TSH))
anTSH<-aov(sqrt(TSH)~Group)
summary(anTSH)

library(MASS)
TukeyHSD(anTSH)