
#loading library
library(Rcpp);library(RcppArmadillo);library(coda);library(lattice);

# BVSG Output

BVSGR <-function(X, Y,Y0, a_init, b_init, r_init, Beta_init, Lambda_init, sigma2_init,ahyper,bhyper,rhyper, Nburn, Niter,Nchain,Nthin){

sourceCpp("BVSG_binary.cpp");

BVSGR<-vector("list",Nchain)
names(BVSGR) <-paste('MCMC_Chain',seq(1:Nchain))

for (i in 1: Nchain) {
BVSGR[[i]] = BVSG(X, Y,Y0, a_init, b_init, r_init, Beta_init, Lambda_init, sigma2_init,ahyper,bhyper,rhyper, Nburn, Niter,Nthin);
names(BVSGR[[i]])=c("Betahat","Lambdahat","ahist","bhist","rhist","sigma2hist","Betahist");
}
return(BVSGR);



}
# BVSG Result Diagnostic

BVSG_Diagnostic<-function(BVSG,Nburn,Niter,Nchain){

		Diagnost <-BVSG[[1]]
		a_mcmc <-Diagnost$ahist ;
		b_mcmc <-Diagnost$bhist ;
		r_mcmc <-Diagnost$rhist ;
		sigma2_mcmc <-Diagnost$sigma2hist ;
		beta_mcmc <-Diagnost$Betahist;


		Diagnost <-BVSG[[2]] ;

		a1_mcmc <-Diagnost$ahist ;
		b1_mcmc <-Diagnost$bhist ;
		r1_mcmc <-Diagnost$rhist ;
		sigma21_mcmc <-Diagnost$sigma2hist ;
		beta1_mcmc <-Diagnost$Betahist;


		a_mcmc <-mcmc.list(mcmc(a_mcmc),mcmc(a1_mcmc));
		b_mcmc <-mcmc.list(mcmc(b_mcmc),mcmc(b1_mcmc));
		r_mcmc <-mcmc.list(mcmc(r_mcmc),mcmc(r1_mcmc));
		sigma2_mcmc <-mcmc.list(mcmc(sigma2_mcmc),mcmc(sigma21_mcmc));
		beta_mcmc <-mcmc.list(mcmc(beta_mcmc),mcmc(beta1_mcmc));
	
BetaPosteriorMean <- matrix(0,dim(Diagnost$Betahist)[2], Nchain)
BetaPosteriorMedian <- matrix(0,dim(Diagnost$Betahist)[2], Nchain)

for ( i in 1: Nchain) {
BetaPosteriorMean[,i] = apply(BVSG[[i]]$Betahist,2,mean)
names(BetaPosteriorMean[,i]) = paste( "BetaPosteriorMean_Chain",i)


BetaPosteriorMedian[,i] = apply(BVSG[[i]]$Betahist,2,median)
names(BetaPosteriorMedian[,i]) = paste( "BetaPosteriorMedian_Chain",i)
}

Cred_Interval = HPDinterval(beta_mcmc,prob=0.95) ;
names(Cred_Interval) = "Credible Interval" ;

split.screen(figs=c(2,1));
split.screen(figs=c(1,3),screen=1);
split.screen(figs=c(1,2),screen=2);

screen(6)
gelman.plot(beta_mcmc,confidence=0.95,autoburnin=FALSE,auto.layout=FALSE,ask=FALSE,xlab="Iterations",ylab="Shrink factor of beta");
screen(3)
gelman.plot(a_mcmc,confidence=0.95,autoburnin=FALSE,auto.layout=FALSE,ask=FALSE,xlab="Iterations",ylab="Shrink factor of a");
screen(4)
gelman.plot(b_mcmc,confidence=0.95,autoburnin=FALSE,auto.layout=FALSE,ask=FALSE,xlab="Iterations",ylab="Shrink factor of b");
screen(5)
gelman.plot(r_mcmc,confidence=0.95,autoburnin=FALSE,auto.layout=FALSE,ask=FALSE,xlab="Iterations",ylab="Shrink factor of r");
screen(7)
gelman.plot(sigma2_mcmc,confidence=0.95,autoburnin=FALSE,auto.layout=FALSE,ask=FALSE,xlab="Iterations",ylab="Shrink factor of sigma2");

Pooled_STD <-summary(sigma2_mcmc)[[1]]
ResultDiagnost<-list(gelman.diag(a_mcmc),gelman.diag(b_mcmc),gelman.diag(r_mcmc),gelman.diag(sigma2_mcmc),gelman.diag(beta_mcmc),BetaPosteriorMean,BetaPosteriorMedian,Pooled_STD,Cred_Interval)
names(ResultDiagnost) <-c( "Potential Scale Reduction Factor for a ", "Potential Scale Reduction Factor for b ", "Potential Scale Reduction Factor for r ","Potential Scale Reduction Factor for sigma2 ","Potential Scale Reduction Factor for Beta ", "Posterior Mean of Beta","Posterior Median of Beta","Pooled Standard Deviance","Confidence Interval of Betahat")
return (ResultDiagnost);

}

myheatmap <- function(X,Lambda) {

### Plotting the Heat Map of Variance Covariance Matrix


label_vector <-c(seq(1:dim(Lambda)[1] ))

mymat      <-solve(t(as.matrix(X) )%*%as.matrix(X) + Lambda)
mysds      <-sqrt(diag(mymat))
mycorrmat  <-mymat/outer(mysds,mysds,FUN="*")

source("myImagePlot.r")

myImagePlot(abs(mycorrmat))


gplot(abs(mycorrmat),edge.lwd= abs(mycorrmat*20),gmode="graph",pad=0.3,loop.cex=5,thresh=0.1,displaylabels=TRUE,label=label_vector,label.bg ="gray90",edge.col="blue",xlab="Dependence Structure among Covariates")


}

CTable <- function( X,Y0, Beta) 
{
P =X%*%Beta ;
P = pnorm(P) ;
n=length(P) ;
Out = matrix(0,n,1) ;
Ctable=matrix(0,n,3);
for( i in 1:length(n) ) {
Out[i,1]=rbinom(1,1,P[i]) ;
}

CTable= cbind(Y0,Out[,1],rep(1,n)) ;
tp =sum(Ctable[CTable[,1]==1 & CTable[,2]==1,3])
fp =sum(Ctable[CTable[,1]==1 & CTable[,2]==0,3])
tn =sum(Ctable[CTable[,1]==0 & CTable[,2]==0,3])
fn =sum(Ctable[CTable[,1]==0 & CTable[,2]==1,3])

ctable <-c(tp,fp,tn,fn) ;
return (ctable);
}