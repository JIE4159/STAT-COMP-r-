f<-file.choose()
data <- read.csv(f,header=TRUE)
z1=data[1:4,3:6]
#####################question  A #########################
## Defining the Gaussian kernel
rbf_kernel <- function(x1,x2,gamma){
  K<-exp(-(1/gamma^2)*t(x1-x2)%*%(x1-x2))
  return(K)
}

###compute H matrix by h=k+cI
kerm<-function(x,gamma,c){
  N=nrow(x)
  km<-matrix(0,N,N)
  Hm<-matrix(0,N,N)
  x<-as.matrix(x)
  for(i in 1:N){
    for(j in 1:N){
      km[i,j]<-rbf_kernel(x[i,],x[j,],gamma)
    }
  }
  Hm=km+c*diag(N)
  return (Hm)
}

H_44<-kerm(z1,4,4)
#####find k vector
kvector<-function(xnew,x,gamma){
  N=nrow(x)
  kvector<-matrix(0,N,1)
  for (i in 1:N){
    kvector[i]<-exp(-(1/gamma^2)*as.matrix((xnew-x[i,]))%*%t(xnew-x[i,]))
  }
  return (kvector)
}
xnew<-data[5,3:6]
k4<-kvector(xnew,z1,4)

####define H-inverse matrix after updating 
H_inv<-function(H,k,c){
  knew<-k
  anew<-solve(H)%*%knew
  rnew<-1+c-t(knew)%*%anew
  Hinv_new1<-rnew[1,1]*(solve(H))+anew%*%t(anew)
  upper<-cbind(Hinv_new1,-anew)
  lower<-cbind(t(-anew),1)
  H_inver_update<-rbind(upper,lower)
  return ((1/rnew[1,1])*H_inver_update)
}
H_inv4=H_inv(H_44,k4,4)

####################  question b #################################
alpha<-function(Hinv){
  N<-nrow(Hinv)
  e<-matrix(1,N,1)
  part1<-2-t(e)%*%Hinv%*%e
  part2<-t(e)%*%Hinv%*%e
  alpha<-(1/2)*Hinv%*%(e+(part1/part2)[1,1]*e)
  return (alpha)
}
alpha5<-alpha(H_inv4)[5]
################ question c #####################################
###add sixth obs 
z2=data[1:5,3:6]
H_55<-kerm(z2,4,4)
xnew6<-data[6,3:6]
k5<-kvector(xnew6,z2,4)
H_inv55=H_inv(H_55,k5,4)
alpha(H_inv55)[6]
####add seventh obs
z3=data[1:6,3:6]
H_66<-kerm(z3,4,4)
xnew7<-data[7,3:6]
k6<-kvector(xnew7,z3,4)
H_inv66=H_inv(H_66,k6,4)
alpha(H_inv66)[7]
#####actual alpha 
zz=data[1:5,3:6]
H55<-kerm(zz,4,4)
H_inv_55<-solve(H55)
alpha(H_inv_55)[5]

zz=data[1:6,3:6]
H66<-kerm(zz,4,4)
H_inv_66<-solve(H66)
alpha(H_inv_66)[6]

zz=data[1:7,3:6]
H77<-kerm(zz,4,4)
H_inv_77<-solve(H77)
alpha(H_inv_77)[7]


############### question d ###################
##### R square three part
KF3<-function(a,x,gamma){
  N<-nrow(x)
  x<-as.matrix(x)
  kf3<-matrix(0,N,N)
  for (i in 1:N){
    for (j in 1:N){
      kf3[i,j]<-a[j]*a[i]*rbf_kernel(x[i,],x[j,],gamma)
    }
  }
  return (sum(kf3))
}

KF33<-KF3(alpha5,z2,4)

KF2<-function(a,x,gamma){
  N<-nrow(x)
  x<-as.matrix(x)
  kf2<-matrix(0,N,N)
  for (i in 1:N){
    for (j in 1:N){
      kf2[i,j]=-2*a[j]*rbf_kernel(x[i,],x[j,],gamma)
    }
  }
  return (sum(kf2)/N)
}
KF22<-KF2(alpha5,z2,4)

Rsquare<-function(kff,kfff){
  rsquare<-kff+kfff+1
  return (rsquare)
}

Rsquare5<-Rsquare(KF33,KF22)

####  distance
dis<-function(x,xtrain,a,gamma){
  for (i in 1:5){
      dis<-a[i]%*%exp(-(1/gamma^2)*as.matrix((x-xtrain[i,]))%*%t(x-xtrain[i,]))
  return (sum(dis))
  }
}

Rsquare5<-Rsquare(KF33,KF22)
dis5<-dis(xnew6,z2,alpha5,4)
dis6<-dis(xnew7,z2,alpha(H_inv55),4)
