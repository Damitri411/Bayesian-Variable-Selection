#Library
library(rjags)
library(R2jags)
library(jagsUI)
library(MASS)
library(base)
library(stats)
library(dplyr)
library(data.table)
library(survival)
library(splines)
library(lme4)
library(lqmm)
path="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS/Horse_shoe"
#Total id
n=5000
X=data.frame(id=1:n,
             X1=runif(n,1,2),
             X2=rnorm(n,-2,1),
             X3=sample(10,size=n, replace=T)/10,
             X4=rexp(n,1),
             X5=rbinom(n,1,0.2)
             )
X$X6=ifelse(X$X1<1.5,rnorm(1,4,3),rnorm(1,-2,1))
X$X7=X$X2^3+X$X4^2
X$X8=sapply(1:nrow(X), function(u) rpois(1,X$X1[u]))

#Beta
beta=c(-1,1,0.02,0,1.5,3,-0.4,0)
beta_0=0.5

X$y=sapply(1:nrow(X), function(u) as.numeric(beta_0+ (beta)%*%t(X[u,2:9])) ) +rnorm(n,0,3)

y=X$y
sink(paste(path,"HS_reg.model.txt", sep="/"))
cat("
model {

for(i in 1:n) 
{
 
 #Regression
 #---------------------
 mu[i]= beta_0 + inprod(beta,X[i,2:9])
 y[i] ~dnorm(mu[i],tau)

}#loop of i

#Prior
#-------------------------------------
  #Error Variance
  
  tau ~ dgamma(0.001,0.001)
  sig = 1/tau
  
  
  #Beta
  #-----
  beta_0 ~dnorm(0,0.001)
  for(p in 1:n_P)
  {
  prec[p]=1/((nu[p]*eta)*(nu[p]*eta))
  beta[p] ~dnorm(0,prec[p])
  }
  
  
  #Horse Shoe Hyper Prior 
  for(p in 1:n_P)
  {
  nu[p] ~ dt(0,1,1)T(0,)
  }
  eta ~dt(0,1,1)T(0,)
  
 
  
  
}
",fill = TRUE)
sink()

sink(paste(path,"HS_reg1.model.txt", sep="/"))
cat("
model {

for(i in 1:n) 
{
 
 #Regression
 #---------------------
 mu[i]= beta_0 + inprod(beta,X[i,2:9])
 y[i] ~dnorm(mu[i],tau)

}#loop of i

#Prior
#-------------------------------------
  #Error Variance
  
  tau ~ dgamma(0.001,0.001)
  sig = 1/tau
  
  
  #Beta
  #-----
  beta_0 ~dnorm(0,0.001)
  for(p in 1:n_P)
  {
  prec[p]=1/(nu1[p]*eta1)
  beta[p] ~dnorm(0,prec[p])
  }
  
  
  #Horse Shoe Hyper Prior 
  for(p in 1:n_P)
  {
  nu1[p] ~ df(1,1)
  nu[p] = pow(nu1[p],0.5)
  
  }
  eta1 ~ df(1,1)
  eta = pow(eta1 , 0.5)
  
 
  
  
}
",fill = TRUE)
sink()


#Data list for Lagged JM
dat_list_HS_reg =list('n'=n,'X'=X , 'y'=y , 'n_P'=ncol(X[,2:9]))

Params_save=c('beta_0', 'beta' ,'sig' , 'eta' ,'nu')



#Run Lagged JM
jagsfit_HS_reg<- jags(data=dat_list_HS_reg, 
                          inits=NULL, 
                          parameters.to.save=Params_save,
                          model.file=paste(path,'HS_reg.model.txt', sep="/"),parallel = TRUE,
                          n.chains=2, n.iter=5000, n.burnin=1000,
                          n.thin=5,n.adapt = 200,
                          DIC=TRUE)


jagsfit_HS_reg1<- jags(data=dat_list_HS_reg, 
                      inits=NULL, 
                      parameters.to.save=Params_save,
                      model.file=paste(path,'HS_reg1.model.txt', sep="/"),parallel = TRUE,
                      n.chains=2, n.iter=5000, n.burnin=1000,
                      n.thin=5,n.adapt = 200,
                      DIC=TRUE)





hs=data.frame(beta_names=paste("beta_",c(1:length(beta)), sep=""),
              Actual_beta=beta,
              True_Importance=ifelse(abs(beta)==0,"No","Yes"),
              HS_importance=jagsfit_HS_reg$mean$nu,
              Estimated_beta=jagsfit_HS_reg$mean$beta,
              LCI=jagsfit_HS_reg$q2.5$beta,
              UCI=jagsfit_HS_reg$q97.5$beta)


hs1=data.frame(beta_names=paste("beta_",c(1:length(beta)), sep=""),
              Actual_beta=beta,
              True_Importance=ifelse(abs(beta)==0,"No","Yes"),
              HS_importance=jagsfit_HS_reg1$mean$nu,
              Estimated_beta=jagsfit_HS_reg1$mean$beta,
              LCI=jagsfit_HS_reg1$q2.5$beta,
              UCI=jagsfit_HS_reg1$q97.5$beta)

library(ggplot2)
library(gridExtra)
p1<-ggplot(data=hs, aes(x=beta_names, y=Estimated_beta))+
 geom_point(size=2, col="red")+
  geom_point( aes(x = beta_names, y = Actual_beta) , shape=1, color="blue", size=2)+
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept=0)+
  theme_bw()

p2<-ggplot(data=hs, aes(x=beta_names, y=HS_importance, group=True_Importance))+
  geom_point(size=3, aes(col=True_Importance))+
  geom_hline(yintercept=0)+
  theme_bw()

p12<-grid.arrange(p1, p2, nrow = 1,
             top = "HS using Cauchy")


p3<-ggplot(data=hs1, aes(x=beta_names, y=Estimated_beta))+
  geom_point(size=2, col="red")+
  geom_point( aes(x = beta_names, y = Actual_beta) , shape=1, color="blue", size=2)+
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept=0)+
  theme_bw()

p4<-ggplot(data=hs1, aes(x=beta_names, y=HS_importance, group=True_Importance))+
  geom_point(size=3, aes(col=True_Importance))+
  geom_hline(yintercept=0)+
  theme_bw()

p34<-grid.arrange(p3, p4, nrow = 1, top="HS using F-dist")

grid.arrange(p12,p34, nrow=2)
