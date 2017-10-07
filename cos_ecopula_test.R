# test 2 different ways to create an empirical copularm(list = ls())
cat("\014")

library(copula)


x=runif(500,-1,1)
y=sin(2*3.14159*x)

################### METHOD 1 ######################################
n=length(x)
vect <- order(x , y )
x1=x
y1=y
for ( i in 1:n) {
  x1[i]=x[vect[i]]
  y1[i]=y[vect[i]]
}
x=x1
y=y1
data<-list(u=x,v=y)
pseudodata <- sapply(data, rank, ties.method = "random") / (n+1 )
u=pseudodata[,1]
v=pseudodata[,2]
cop=u
for(i in 1:n) {
  cpt=0
  for(j in 1:i) {
    if (( u[i]>=u[j]) && ( v[i]>=v[j])) {
      cpt=cpt+1
    }
  }
  cop[i]=cpt/(n + 1)
}

################### METHOD 2 ######################################
XX = cbind(x,y)
ec <- C.n(pseudodata, XX)

par(mfrow=c(1,2)) 
plot(cop)
plot(ec)