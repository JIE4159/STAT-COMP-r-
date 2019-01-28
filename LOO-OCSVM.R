##### import dataset and standardize predictors
f<-file.choose()
data <- read.table(f,header = TRUE)
zvar<-data.matrix(data[1:61,2:5])
Stzvar<-(zvar-mean(zvar))/sd(zvar)
y<-data[1:61,1]
pbdata<-data.frame(y,Stzvar)
##########  set sigma=4 C=50  #########################
####### start a Initial C making zero misclassifcation rate####
#### km matrix##
rbf_kernel <- function(x1,x2,sigma){
  K<-exp(-(1/gamma^2)*t(x1-x2)%*%(x1-x2))
  return(K)
}
km<-function(x,sigma,c){
  N=nrow(x)
  km<-matrix(0,N,N)
  x<-as.matrix(x)
  for(i in 1:N){
    for(j in 1:N){
      km[i,j]<-rbf_kernel(x[i,],x[j,],sigma)
    }
  }
  return (km)
}
################alpha vector##############
alphas<-function(x,sigma,C){
  len=dim(x)[1]+1
  A=matrix(0,len,len)
  A[1:len-1,1:len-1]=km(x,sigma,C)+identity(len-1)/C
  A[1:len-1,len]=matrix(1,len-1,1)
  A[len,1:len-1]=matrix(1,1,len-1)
  A[len,len]=0
  B=matrix(0,len,1)
  B[len][1]=1
  result=solve(A,B)
  result1=data.frame(result,A)
  return (result1)
  }

clas<-function(a,b,x,sigma,C){
  cls<-numeric(30)
  for (i in 1:30){
    for (j in 1:30){
      cls[i]=b+cumsum(a[j]*km(x,sigma,C)[i,j])
    }
  }
  return (cls)
}
normy<-function(cls,C){
  ymax=max(cls)
  ymin=min(cls)
  normy<-1-2*abs(cls)/(C*sign(cls)*(ymax+ymin)+(ymax-ymin)) ##C=50 makes misclassification=0
  return (sign(normy))
}
a1<-alphas(Stzvar[1:30,],4,50)[1:30,1]
b1<-alphas(Stzvar[1:30,],4,50)[31,1]
A<-alphas(Stzvar[1:30,],4,50)[,2:32]
res=clas(a1,b1,Stzvar,4,50)
normy(res,50) ##satisfy zero misclassification
#######grid search while loo error not eqal to 0 ###
y_loo<-function(sigma,C){
  a11<-alphas(Stzvar[1:30,],sigma,C)[1:30,1]
  A1<-alphas(Stzvar[1:30,],sigma,C)[,2:32]
  y_loo1=numeric(30)
  for (i in 1:30){
    y_loo1[i]<- -(a11[i])/(solve(A1))[i,i]
  }
  res<-normy(y_loo1,C)
  return (res)
}

####### find best C######
best.C<-function(sigma){
  com<-rep(1,30)
  C<-50
  while (identical(y_loo(sigma,C),com))
  {
      C=0.8*C
  }
  return (C)
}
best.C(sigma=4)
best.C(sigma=0.1)


################# problem 2 ####################################
#### compute omega, zt, ###########
data2=pbdata[1:30,-1]
y=data.frame(data2)
omega=cov(data2)
omega_inv=solve(omega)
x=chol(omega_inv)
##create dataset for standarized data and ########
###response are z, convariate is x matrix ###
z=matrix(0,4,30)
for (i in 1:30){
  z[,i]=x%*%t(y[i,])
}
###################variables selection by linear regression to choose miu#### 
lmresult<-function(x,y){
  miu_star<-matrix(0,30,5)
  for (i in 1:30){
  z[,i]=x%*%t(y[i,])
  data<-data.frame(z[,i],x) 
  lm=lm(z...i.~.,data = data)
  miu_star[i,]=lm$coefficients
  }
 return (miu_star)
}
miu_star<-lmresult(x,y)

#######################  process monitoring ###############################
C=1
lambda<-numeric(30)
trigger<-numeric(30)
for (i in 1:30){
  yt<-as.matrix(y[i,])
  miut<-as.matrix(miu_star[i,1:4])
  lambda[i]<-2*yt%*%omega_inv%*%miu_star[i,2:5]-t(miu_star[i,2:5])%*%omega_inv%*%miu_star[i,2:5]
  trigger[i]<-lambda[i]-C
  print(c(i,trigger[i]))
}
sum(sign(trigger)==1)
rates<-sum(sign(trigger)==1)/30



#################### problem 3 dual quadratic SVDD ################################
library(quadprog)
Data <- as.matrix(data, ncol=5)
Xtrain <- as.matrix(Data[1:30,2:5])

require('quadprog')
## Defining the Gaussian kernel

svmtrain <- function(X,C=Inf, gamma,esp=1e-10){
  
  N<-nrow(X)
  Dm<-matrix(0,N,N)
  X<-as.matrix(X)
  for(i in 1:N){
    for(j in 1:N){
      Dm[i,j]<-2*rbf_kernel(X[i,],X[j,],gamma)
    }
  }
  Dm<-Dm+diag(N)*1e-12 # adding a very small number to the diag, some trick
  
  dv<-t(rep(1,N))
  meq<-1
  Am<-cbind(matrix(1,N),diag(N))
  bv<-rep(0,1+N) # the 1 is for the sum(alpha)==0, others for each alpha_i >= 0
  Am<-cbind(Am,-1*diag(N))
  bv<-c(cbind(matrix(bv,1),matrix(rep(-C,N),1)))
  alpha_org<-solve.QP(Dm,dv,Am,meq=meq,bvec=bv)$solution
  return (alpha_org)
  
}


alpha<-svmtrain(Xtrain,5, 4,esp=0)
alphaindx<-which(alpha>0,arr.ind=TRUE)
alpha11<-alpha[alphaindx]
nSV<-length(alphaindx)
### Predict the class of an object X
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
kf33<-KF3(alpha,Xtrain,4)

KF2<-function(a,x,gamma){
  N2<-length(alphaindx)
  N<-nrow(x)
  x<-as.matrix(x)
  kf2<-matrix(0,N2,N)
  for (i in 1:N2){
    for (j in 1:N){
      kf2[i,j]=-2*a[j]*rbf_kernel(x[alphaindx[i],],x[j,],gamma)
    }
  }
  return (sum(kf2)/N2)
}

kf22<-KF2(alpha,Xtrain,4)

Rsquare<-function(kff,kfff){
  rsquare<-kff+kfff
  return (rsquare)
}

Rsquare5<-Rsquare(kf33,kf22)

##############  distance############

Xtest <- as.matrix(Data[31:61,2:5])


D<-numeric(31)
dis<-numeric(31)
for (i in 1:31){
  for (j in 1:30){
    dis[i]=sum(-2*alpha[j]*rbf_kernel(Xtest[i,],Xtrain[j,],4))
    D[i]=dis[i]+kf33
  }
}

errorates<-function(d,r){
  count<-0
  for (i in 1:31){
    if (d[i]<=r){
      count=count+1
    }
  }
  return (count/31)
}
errorates(D,Rsquare5)

##################### problem 2   ##################
####step 1 try sigma=5, C=o.5 ######
newdata <- as.matrix(Data[21:50,2:5])
y<-data[21:50,1]
class1<-newdata[1:10,]
class2<-newdata[11:30,]
alpha11<-svmtrain(class1,0.5, 4,esp=1e-10)
alpha22<-svmtrain(class2,0.5, 4,esp=1e-10)
first_layer<-c(alpha11,alpha22)
######## step 2 compare radius and confidence measurement#######
kf331<-KF3(alpha11,Xtrain,4)
kf221<-KF2(alpha11,Xtrain,4)
Rsquare1<-Rsquare(kf331,kf221)
Radius1<-Rsquare5+1

kf332<-KF3(alpha22,Xtrain,4)
kf222<-KF2(alpha22,Xtrain,4)
Rsquare2<-Rsquare(kf332,kf222)
Radius2<-Rsquare5+1
##confidence measurement 
##class 1
center1<-c(mean(newdata[1:10,1]),mean(newdata[1:10,2]),mean(newdata[1:10,3]),mean(newdata[1:10,4]))
CM1<-numeric(10)
for (i in 10){
  CM1[i]<-exp(-(sum((center1-class1[i,])^2))/Radius1)
  }
## calss 2
center2<-c(mean(newdata[11:30,1]),mean(newdata[11:30,2]),mean(newdata[11:30,3]),mean(newdata[11:30,4]))
CM2<-numeric(20)
for (i in 20){
  CM2[i]<-exp(-(sum((center1-class2[i,])^2))/Radius1)
} 
