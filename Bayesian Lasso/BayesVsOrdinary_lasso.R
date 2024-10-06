
#Attach libraries
library(care)
library(glmnet)
library(HDCI)
library(monomvn)
library(plotrix)
library(graphics)
library(statmod)
library(MASS)
library(invgamma)
library(base)


#Data formation

data(efron2004)
d_x=efron2004$x
d_y=efron2004$y

colnames(d_x)

dat=data.frame(age=as.numeric(),
               sex=as.numeric(),
               bmi=as.numeric(),
               bp=as.numeric(),
               s1=as.numeric(),
               s2=as.numeric(),
               s3=as.numeric(),
               s4=as.numeric(),
               s5=as.numeric(),
               s6=as.numeric(),
               y=as.numeric())
for(i in 1:nrow(d_x))
{
  dat[i,]=c(d_x[i,],d_y[i])
  i=i+1
}
#Standardize
dat1=dat
for( i in 1:ncol(dat))
{
  dat1[,i]=(dat1[,i]-mean(dat1[,i]))/sd(dat1[,i])
  i=i+1
}
sd_dat1=apply(dat1,2,sd)
mean_dat1=apply(dat1,2,mean)
dat=dat1

#########################################################################################################
#Click the link below
#https://drsimonj.svbtle.com/ridge-regression-with-glmnet
############################################################################################################
#                                              Trace plot for lasso                                        #
############################################################################################################

#nfold cv
cv_lasso=cv.glmnet(x=as.matrix(dat[,1:10]), y=dat$y ,alpha = 1 , nfolds = 20)

#cv_lasso$lambda.min
#cv_lasso$lambda.1se

plot(cv_lasso)
text(x=log(c(cv_lasso$lambda.min,cv_lasso$lambda.1se)) , y=c(0.6,0.6) ,labels=c("l_min","l_1se") )
coef_cv=coef(cv_lasso , s="lambda.min")
bet_cv=t(as.matrix(coef_cv))


lambda_lasso=c(exp(seq(0,-9,-0.1)),0)
#ordinary lasso
fit_olso=glmnet(x=as.matrix(dat[,1:10]), y=dat$y ,family = "gaussian" ,lambda=lambda_lasso ,
                alpha = 1 ,standardize = T , intercept = T)


#betas
#fit_olso$beta

#intercept
#fit_olso$a0



bet=NULL
bet=t(as.matrix(coef(fit_olso)))
bet1=abs(bet)
l1_norm=apply(bet1,1,sum)
l1_rel=l1_norm/max(l1_norm)
bet2_olso=cbind(bet,l1_rel)

par(mfrow=c(1,3))

matplot(bet2_olso[,12],bet2_olso[,2:11],type="l",col=3:12,lty=1,xlab="Relative L1",ylab="Standardized Coefficients",main="Lasso" )
text(x=rep(bet2_olso[nrow(bet2_olso),12],10),y=bet2_olso[nrow(bet2_olso),2:11]+0.005,labels=colnames(bet2_olso)[2:11] , col=3:12)

l1_cv_rel=sum(abs(bet_cv[2:11]))/max(l1_norm)
abline(v=l1_cv_rel , lty=2)
abline(h=0 , lty=2)
text(x=l1_cv_rel-0.02 , y=-0.41 , labels=c("cv"))

#trace plot for lasso
#plot(fit_olso,xvar="lambda",label=T)

#sequence of lambdas
#lambda_lass=fit_olso$lambda


############################################################################################################
#                                              Trace plot for blasso                                        #
############################################################################################################
#sd_dat=apply(dat,2,sd)
#mean_dat=apply(dat,2,mean)


pos_est=NULL
lambda_blass=c(exp(seq(7,-7,-0.1)))
for(j in 1:length(lambda_blass))
{
  fit_bl=blasso(X=as.matrix(dat[,1:10]), y=dat$y, T = 1000, thin = 100, RJ = FALSE, M = NULL,
                    beta = NULL, lambda2 = lambda_blass[j], case = c("default"), s2 = 1,
                    rd = FALSE, ab = NULL, rao.s2 = TRUE,
                    normalize = TRUE)
  pos_med=c(as.vector(apply(fit_bl$beta,2,median)),as.vector(apply(fit_bl$tau2i,2,median)),as.numeric(median(fit_bl$s2)),lambda_blass[j])
  pos_est=rbind(pos_est,pos_med)
  
  j=j+1
  
}
pos_est1=as.matrix(pos_est)
colnames(pos_est1)=c(colnames(dat1)[1:10],paste0("tau_",c(1:10),sep=""),"sig_sq","lambda")
rownames(pos_est1)=c(1:nrow(pos_est1))

bet_bl=NULL
bet_bl=pos_est1[,c(1:10,22)]
bet_bl=rbind(bet_bl,c(as.vector(coef(lm(y~.,data=dat1))[2:11]),0))
bet1_bl=abs(bet_bl)
l1_bl_norm=apply(bet1_bl[,1:10],1,sum)
l1_rel=l1_bl_norm/max(l1_bl_norm)
bet2_bl=cbind(bet_bl,l1_rel)



matplot(bet2_bl[,12],bet2_bl[,1:10],type="l",col=3:12,lty=1,xlab="Relative L1",ylab="Standardized Coefficients",main="Bayesian Lasso" )
text(x=rep(bet2_bl[nrow(bet2_bl),12],10),y=bet2_bl[nrow(bet2_bl),1:10]+0.009,labels=colnames(bet2_bl)[1:10] , col=3:12)
abline(h=0 , lty=2)


#bayesian lasso along with sampling of lambda

#MML
X=as.matrix(dat1[,1:10])
y=as.matrix(dat1$y)

fit_reg=lm(y~. , data=dat1)
coeff_reg=coef(fit_reg)
sig2_reg=sum((fit_reg$residuals)^2)/(fit_reg$df.residual)

l_0=sqrt(sig2_reg)/mean(abs(coeff_reg[2:11]))
#l_0=0.3
tau2_0=rexp(10, rate=l_0^2/2)
sig2_0=sig2_reg
b_0=mvrnorm(1,rep(0,10),sig2_0*diag(tau2_0))


post_est=t(as.matrix(c(b_0,tau2_0,sig2_0,l_0)))
colnames(post_est)=c(colnames(dat1)[1:10], paste("tau2_",1:10,sep=""),"sig_sq","lambda")
post_est=data.frame(post_est)

iter_max=1000
burn_in=100
for(i in 2:(iter_max+1))
{
  b=post_est[i-1,1:10]
  t2=post_est[i-1,11:20]
  s2=post_est[i-1,21]
  l=post_est[i-1,22]
  
  
  b=mvrnorm(1, solve(t(X)%*%X +diag(t2))%*%t(X)%*%y ,s2*solve(t(X)%*%X +diag(t2)))
  
  s2=rinvgamma(1,shape=(442-1+10)/2 , scale=1)*(as.numeric(t(y-X%*%b)%*%(y-X%*%b)+ t(b)%*%diag(t2)%*%b))/2
  
  t2=sapply(1:10 , function(x) 1/rinvgauss(1,sqrt(l^2*s2/b[x]^2),l^2))
  
  fit_int=blasso(X=as.matrix(dat[,1:10]), y=dat$y, T = 1000, thin = 100, RJ = FALSE, M = NULL,
                 beta = b, lambda2 = (l^2), case = c("default"), s2 = s2,
                 rd = FALSE, ab = NULL, rao.s2 = FALSE,
                 normalize = TRUE)
  
  l=sqrt((2)/mean(as.vector(apply(1/fit_int$tau2i,2,mean))))
  #l=sqrt(rgamma(1,shape=11,scale=1/(sum(t2)/2 )))
  
  
  post_est[i,]=c(b,t2,s2,l)
  i=i+1
  
}

post_med=apply(post_est[(burn_in+2):nrow(post_est),],2,median)

#Lambda by MML
lam_med=post_med[22]
l1_mml_rel=sum(abs(post_med[1:10]))/max(l1_bl_norm)
abline(v=l1_mml_rel , lty=2)
text(x=l1_mml_rel-0.02 , y=-0.41 , labels=c("mml"))


delta=1/(10*(lam_med)^2)
#With prior
fit_blasso=blasso(X=as.matrix(dat[,1:10]), y=dat$y, T = 1000, thin = 100, RJ = FALSE, M = NULL,
               beta = NULL, lambda2 = 0.1, case = c("default"), s2 = 1,
               rd = c(1,delta), ab = NULL, rao.s2 = TRUE,
               normalize = TRUE)
lam_med_prior=median(sqrt(fit_blasso$lambda2))

l1_bl_pr_rel=sum(abs(apply(fit_blasso$beta,2,median)))/max(l1_bl_norm)
abline(v=l1_bl_pr_rel , lty=2, col="red")
text(x=l1_bl_pr_rel-0.03 , y=-0.3 , labels=c("pr"), col="red")

############################################################################################################
#                                              Trace plot for ridge                                        #
############################################################################################################
lambda_ridge=c(exp(seq(9,-4,-0.1)),0)

#ordinary ridge
fit_oridge=glmnet(x=as.matrix(dat[,1:10]), y=dat$y ,family = "gaussian" ,lambda=lambda_ridge, 
                  alpha = 0 ,standardize = T , intercept = T)

bet_r=NULL
bet_r=t(as.matrix(coef(fit_oridge)))
bet_r1=abs(bet_r)
l1_norm_r=apply(bet_r1,1,sum)
l1_rel_r=l1_norm_r/max(l1_norm_r)
bet_r2=cbind(bet_r,l1_rel_r)



matplot(bet_r2[,12],bet_r2[,2:11],type="l",col=3:12,lty=1,xlab="Relative L1",ylab="Standardized Coefficients",main="Ridge" )
text(x=rep(bet_r2[nrow(bet_r2),12],10),y=bet_r2[nrow(bet_r2),2:11]+0.005,labels=colnames(bet_r2)[2:11] , col=3:12)
abline(h=0 , lty=2)




############################################################################################################
#                                           CI and estimate plots code                                    #
############################################################################################################

#y1=posterior means
#liw=y1-5%ci
#uiw=95%ci-y1
#xlim=c(-max(|all coeff|),max(|all coeff|))




y1=post_med[1:10]
liw=y1-as.vector(sapply(1:10,function(x) quantile(post_est[,x],probs=c(0.025))))
uiw=as.vector(sapply(1:10,function(x) quantile(post_est[,x],probs=c(0.975))))-y1

#posterior estimates of blasso
plotCI(x=y1,y=1:10,uiw=uiw,liw=liw,pt.bg=par("bg"),pch=10,err="x",xlim=c(-1,1),ylim=c(0,11),xlab="Standardized Coefficients", ylab="Variable number",main="CI and estimates")

#In place of x in points(..,x= ,,) attach other estimates
#Regression
points(y=1:10,x=coeff_reg[2:11], col="blue", pch=4)
#n-fold cv lasso n=20
points(y=(1:10)-0.3,x=bet_cv[2:11], col="red", pch=2)
#blasso l1=lasso l1
points(y=(1:10)+0.3,x=as.vector(bet2_olso[as.numeric(which.min(abs(bet2_olso[,12]-l1_mml_rel))),2:11])
, col="green", pch=6)
abline(v=0 , lty=3)

text(x=rep(0.8,10),y=1:10,labels=colnames(post_est)[1:10] , col="black")

############################################################################################################
#                                               Stats Reqd                                                 #
############################################################################################################
#L1 relative to regression for ordinary lasso using cv
l1_cv_rel
#L1 relative to regression for Bayesian lasso using MML
l1_mml_rel
#L1 relative to regression for Bayesian lasso using prior
l1_bl_pr_rel

#Lambda for olso
cv_lasso$lambda.min
cv_lasso$lambda.1se
#Lambda for blasso (estimate and CI) by MML
quantile(post_est$lambda,probs=c(0.025,0.5,0.975))
#Lambda for blasso (estimate and CI) by prior
quantile(sqrt(fit_blasso$lambda2),probs=c(0.025,0.5,0.975))
#prior used
rd=c(1,delta)

y1=c(post_med[1:10],apply(fit_blasso$beta,2,median))
liw=y1-c(as.vector(sapply(1:10,function(x) quantile(post_est[,x],probs=c(0.025)))) ,as.vector(sapply(1:10,function(x) quantile(fit_blasso$beta[,x],probs=c(0.025)))))
uiw=c(as.vector(sapply(1:10,function(x) quantile(post_est[,x],probs=c(0.975)))),as.vector(sapply(1:10,function(x) quantile(fit_blasso$beta[,x],probs=c(0.975)))))-y1

#posterior estimates of blasso
plotCI(x=y1,y=as.vector((c(2*1:10-0.25,2*1:10+0.25))),uiw=uiw,liw=liw,pt.bg=par("bg"),pch=10,col=c(rep(c("red"),10),rep(c("blue"),10)),err="x",xlim=c(-1,1),ylim=c(0,22),xlab="Standardized Coefficients", ylab="Variable number",main="prior Vs MML")
text(x=rep(0.8,10),y=2*1:10,labels=colnames(post_est)[1:10] , col="black")
abline(v=0, lty=2)




#red mml
#blue prior

############################################################################################################
#                                                  End                                                     #
############################################################################################################

