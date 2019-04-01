install.packages("BGLR")
library(BGLR)
data(wheat)
dim(wheat.X) ## 599 individuals and 1279 markers
wheat.X[1:10,1:10] ## use the dominant marker
dim(wheat.Y)
head(wheat.Y) ## four different envrionments, 599 samples. phenotype data
Y<-wheat.Y[,4] ## use the 4th environment data
x<-scale(wheat.X) ## centering data with mean 0 and variation 1
mean(x[,1])
mean(x[,2]) ## after the sacle, the mean of diffenet marker close to 1
var(x[,1])
x[1:10,1:10]

regression<-lm (Y~x[,1]) ## 0.719 is larger than 0.05, not significant marker
summary(regression)

anova(regression) 
anova(regression)[1,5] ## extract the p-value
p<-vector()
for(i in 1:ncol(x)){
  regression<-lm (Y~x[,i])
  p[i]<-anova(regression)[,5]
  
}
head(p)
length(p)
bonferroni<-0.05/ncol(x)
bonferroni
plot(-log(p,base=10),pch=12,cex=0.6,
     ylab=expression(paste("-log"[10],"(p)", sep="")), 
     xlab="DarTs Markers",main= "Gwas_environment 4",ylim=c(0,15))
abline(h=-log(bonferroni,base=10),col="red", lty=2)

######################################################################
D=dist(x, method="euclidean", diag=T, upper=T)
as.matrix(D)
Diss=1-as.matrix(D)
Diss
EVD=eigen(Diss)
plot(EVD$vectors[,1:2])
set.seed(1) ## get the same result everytime
new_df= cbind(EVD$vectors[,1],EVD$vectors[,2])
fit=kmeans(new_df,centers=2)
color=c("blue","red")
plot(EVD$vectors[,1:2], col=color[fit$cluster])

##fit=kmeans(new_df,centers=3) ## change the number of kmeans
##color=c("blue","red","green")
##plot(EVD$vectors[,1:2], col=color[fit$cluster])
#######################################################################

p<-vector()
pop.cofactor=factor(fit$cluster)
pop.cofactor
for(i in 1:ncol(x)){
  regression<-lm (Y~pop.cofactor+x[,i])
  p[i]<-anova(regression)[,5]
  
}
head(p)
length(p)
bonferroni<-0.05/ncol(x)
bonferroni
plot(-log(p,base=10),pch=12,cex=0.6,
     ylab=expression(paste("-log"[10],"(p)", sep="")), 
     xlab="DarTs Markers",main= "Gwas_environment 4",ylim=c(0,15))
abline(h=-log(bonferroni,base=10),col="red", lty=2)
which(p<= bonferroni)















