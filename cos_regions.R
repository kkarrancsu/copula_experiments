rm(list = ls())
cat("\014")

cosdv_regions_2d <- function(x, y) {
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
  ss=cop[1]
  
  regionIdxVec <- vector()
  regionIdx = 1
  
  for ( i in 2:n) {
    ss=c(ss,cop[i])
    if (( is.unsorted (ss))&&( is.unsorted (rev(ss)))) {
      regionIdxVec[regionIdx] = i; 
      regionIdx = regionIdx + 1;
      
      ss=cop[i]
    }
    else {
      if (i==n) {
        regionIdxVec[regionIdx] = i; 
      }
    }
  }
  
  print(length(regionIdxVec))
  print(regionIdxVec)
  
  return(cop)
}

cos_regions_3d <- function(x, y,z) {
  
  n=length(x)
  vect <- order(x,y,z )
  x1=x
  y1=y
  z1=z
  for ( i in 1:n) {
    x1[i]=x[vect[i]]
    y1[i]=y[vect[i]]
    z1[i]=z[vect[i]]
  }
  x=x1
  y=y1
  z=z1
  
  data<-list(u=x,v=y,w=z)
  pseudodata <- sapply(data, rank, ties.method = "random") / (n+1 )
  u=pseudodata[,1]
  v=pseudodata[,2]
  w=pseudodata[,3]
  
  cop=u
  ss=cop[1]
  
  regionIdxVec <- vector()
  regionIdx = 1
  
  for(i in 1:n) {
    cpt=0
    for(j in 1:i) {
      if (( u[i]>=u[j]) && ( v[i]>=v[j]) && ( w[i]>=w[j])) {
        cpt=cpt+1
      }
    }
    cop[i]=cpt/(n + 1)
  }
  
  for ( i in 2:n) {
    ss=c(ss,cop[i])
    if (( is.unsorted (ss))&&( is.unsorted (rev(ss)))) {
      
      # print("new region IN BETWEEN")
      regionIdxVec[regionIdx] = i; 
      regionIdx = regionIdx + 1;
      
      ss=cop[i]
    }
    else {
      if (i==n) {
        # print("new region at N!")
        regionIdxVec[regionIdx] = i;
      }
    }
  }
  
  print(length(regionIdxVec))
  return(cop)
}

n <- 2500
x = runif(n,-4,4)
y = runif(n,-4,4)
# z = y*sin(x) - x*cos(y)
z = -x^2-y^2+6
ecop_3d <- cos_regions_3d(x,y,z)

x2 = runif(500)
y2 = 4*(x2-0.5)^2
ecop_2d <- cosdv_regions_2d(x2,y2)
# plot(ecop_2d)

par(mfrow=c(1,2))
plot(ecop_2d)
plot(ecop_3d)

# plot3D::scatter3D(x,y,z)
# plot(ecop)
# scatterplot3d::scatterplot3d(x,y,z)

# library("car")
# scatter3d(x,y,z)