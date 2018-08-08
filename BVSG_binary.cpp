#include<RcppArmadillo.h>
#include<cstdio>
#include<cstdlib>
#include<iostream>

//Sampling Scheme from Inverse Gaussian Distribution
//[[Rcpp::export]]
Rcpp::NumericVector inversegauss(int k,Rcpp::NumericVector mu_ig_v,Rcpp::NumericVector lambda_ig_v){

//variable declaration
Rcpp::NumericVector x_ig_v(k);
Rcpp::NumericVector UniRand   = Rcpp::runif(k); //sampling a vector from U(0,1)
Rcpp::NumericVector NormRand  = Rcpp::rnorm(k); //sampling a vector from N(0,1)

for(int i =0; i<k ; i++){

    //variable declaration
    double mu_ig = mu_ig_v(i);
    double lambda_ig = lambda_ig_v(i);
    double y_ig;
    double x_ig;

        // returns 0, when we sample from inverse gauss(0.0,lambda)
        if(mu_ig!=0.0) {

            y_ig      = pow(NormRand[i],2) ;
            x_ig      = mu_ig + (pow(mu_ig,2)*y_ig)/(2*lambda_ig) - (mu_ig/(2*lambda_ig))*sqrt(4*mu_ig*lambda_ig*y_ig+ pow(mu_ig,2)*pow(y_ig,2));

                if (UniRand[i] >= mu_ig/(mu_ig+x_ig)) {
                        x_ig_v[i] = (pow(mu_ig,2))/x_ig;
                } else {
                            x_ig_v[i] = x_ig;
                        }

        } else {
            x_ig_v[i] = 0.0;
        }

}
return x_ig_v;
}


//Sampling Beta from Multivariate Normal
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::NumericVector sampleBetainC(double sigma2, Rcpp::NumericMatrix Lambda,Rcpp::NumericMatrix XtX,Rcpp::NumericVector XtYminusYbar){

int n=XtX.nrow(), p=XtX.ncol();
arma::mat xtx(XtX.begin(),n,p,false) ;
arma::mat lambda(Lambda.begin(),p,p,false) ;
arma::vec xtyminusybar(XtYminusYbar.begin(),XtYminusYbar.size(),false);

arma::mat cov_beta  =  arma::inv(xtx + lambda);
arma::vec mean_beta = cov_beta*xtyminusybar;

/*
Rcpp::Environment mvtnorm("package:mvtnorm");
Rcpp::Function rmvnorm =mvtnorm["rmvnorm"];
*/
Rcpp::NumericVector mean = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mean_beta));
Rcpp::NumericMatrix cov  = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(cov_beta));

/*
Rcpp::String Boolean = "svd";
Rcpp::NumericVector beta= rmvnorm(1,Rcpp::Named("mean",mean),Rcpp::Named("sigma",cov),Rcpp::Named("method",Boolean));



Rcpp::Environment myEnv("package::base");
Rcpp::Function  mychol  = myEnv["chol"];
Rcpp::NumericMatrix chol_cov(p,p);
chol_cov = mychol(cov);
*/
arma::vec beta(p);

arma::vec normrand = Rcpp::as<arma::vec>(Rcpp::rnorm(p)) ;
beta = mean_beta + sqrt(sigma2)*arma::chol(cov_beta)*normrand;

return Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(beta));

}


//Sampling from C_matrix
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::NumericMatrix sampleC_matrixinC(Rcpp::NumericVector beta,double sigma2,double r,double a,double b){

//variable declaration
int p=beta.size();
arma::colvec one(p,arma::fill::ones);
arma::colvec Beta(p);
Beta  = Rcpp::as<arma::colvec>(Rcpp::wrap(beta));

arma::mat Beta_Plus(p,p);
Beta_Plus   = Beta*one.t() + one*Beta.t();

arma::mat Beta_minus(p,p);
Beta_minus  = Beta*one.t() - one*Beta.t();

arma::mat Pmat(p,p);
arma::mat onemat(p,p);
onemat.fill(1.0);

Pmat   = r*b*(arma::abs(Beta_Plus) - arma::abs(Beta_minus))*(1/(2*sqrt(sigma2)));
Pmat   = onemat + arma::exp(Pmat);
Pmat   = arma::pow(Pmat,-1.0);

arma::mat unifmat(p,p) ;
unifmat = arma::randu<arma::mat>(p,p);

arma::mat diffmat     = unifmat - Pmat;

arma::mat Cmat        = arma::sign(diffmat);
          Cmat        = arma::symmatu(Cmat)-arma::diagmat(Cmat)+ arma::diagmat(one);


Rcpp::NumericMatrix Cmatrix = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(Cmat));

return Cmatrix ;
}

// norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to be in the interval
// (a,b) via rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double norm_rs(double a, double b)
{
   double  x;
   x = Rf_rnorm(0.0, 1.0);
   while( (x < a) || (x > b) ) x = norm_rand();
   return x;
}

// half_norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) (with a > 0) using half normal rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double half_norm_rs(double a, double b)
{
   double   x;
   x = fabs(norm_rand());
   while( (x<a) || (x>b) ) x = fabs(norm_rand());
   return x;
}

// unif_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using uniform rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double unif_rs(double a, double b)
{
   double xstar, logphixstar, x, logu;

   // Find the argmax (b is always >= 0)
   // This works because we want to sample from N(0,1)
   if(a <= 0.0) xstar = 0.0;
   else xstar = a;
   logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);

   x = R::runif(a, b);
   logu = log(R::runif(0.0, 1.0));
   while( logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
   {
      x = R::runif(a, b);
      logu = log(R::runif(0.0, 1.0));
   }
   return x;
}

// exp_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using exponential rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double exp_rs(double a, double b)
{
  double  z, u, rate;

//  Rprintf("in exp_rs");
  rate = 1/a;
//1/a

   // Generate a proposal on (0, b-a)
   z = R::rexp(rate);
   while(z > (b-a)) z = R::rexp(rate);
   u = R::runif(0.0, 1.0);

   while( log(u) > (-0.5*z*z))
   {
      z = R::rexp(rate);
      while(z > (b-a)) z = R::rexp(rate);
      u = R::runif(0.0,1.0);
   }
   return(z+a);
}


// rnorm_trunc( mu, sigma, lower, upper)
//
// generates one random normal RVs with mean 'mu' and standard
// deviation 'sigma', truncated to the interval (lower,upper), where
// lower can be -Inf and upper can be Inf.
//======================================================================

// [[Rcpp::export]]
double rnorm_trunc (double mu, double sigma, double lower, double upper)
{
int change;
 double a, b;
 double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
 double z, tmp, lograt;

 change = 0;
 a = (lower - mu)/sigma;
 b = (upper - mu)/sigma;

 // First scenario
 if( (a == R_NegInf) || (b == R_PosInf))
   {
     if(a == R_NegInf)
       {
     change = 1;
     a = -b;
     b = R_PosInf;
       }

     // The two possibilities for this scenario
     if(a <= 0.45) z = norm_rs(a, b);
     else z = exp_rs(a, b);
     if(change) z = -z;
   }
 // Second scenario
 else if((a * b) <= 0.0)
   {
     // The two possibilities for this scenario
     if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
       {
     z = norm_rs(a, b);
       }
     else z = unif_rs(a,b);
   }
 // Third scenario
 else
   {
     if(b < 0)
       {
     tmp = b; b = -a; a = -tmp; change = 1;
       }

     lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
     if(lograt <= logt2) z = unif_rs(a,b);
     else if((lograt > logt1) && (a < t3)) z = half_norm_rs(a,b);
     else z = exp_rs(a,b);
     if(change) z = -z;
   }
   double output;
   output = sigma*z + mu;
 return (output);
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::NumericVector sampleTruncY(Rcpp::NumericMatrix X, Rcpp::NumericVector Y0, Rcpp::NumericVector Beta, double sigma2 ){

// variable declarations;
arma::mat x ;
arma::vec beta ;
arma::mat mu ;
int p =Y0.size();
double rand;
double mean ;
Rcpp::NumericVector Z(p);
Rcpp::NumericVector v(p) ;

x    = Rcpp::as<arma::mat>(X) ;
beta = Rcpp::as<arma::vec>(Beta);

mu = x*beta ;
v =Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mu)) ;

for ( int i =0 ; i < p ; i++) {

        mean = v[i] ;
     if ( Y0(i) == 1 )
        {
         while ( rand <= 0 || rand > 100) {
         rand = Rcpp::as<double>(Rcpp::rnorm(1))*sqrt(sigma2) + mean ;
         }
         Z[i] = rand ;
         rand = 0;
        }
        else if(Y0(i) == 0)
        {
         while ( rand  >= 0 || rand <- 100) {
         rand = Rcpp::as<double>(Rcpp::rnorm(1))*sqrt(sigma2) + mean ;
         }
         Z[i] = rand ;
         rand = 0;
        }

            }

return Z ;

}




//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::NumericMatrix sampleLambdainC(Rcpp::NumericVector beta,Rcpp::NumericMatrix C_matrix,double sigma2,double r,double a,double b){

//variable declaration
int p = beta.size();
int d = p*(p-1)/2;
int count ;


Rcpp::NumericMatrix Lambda_off(p,p);
Rcpp::NumericMatrix betam(p,p);
Rcpp::NumericVector b_mu_ig(d);
Rcpp::NumericVector b_lambda_ig(d);
Rcpp::NumericVector b_ig(d);
Rcpp::NumericVector a_mu_ig(p);
Rcpp::NumericVector a_lambda_ig(p);
Rcpp::NumericVector a_ig(p);


arma::mat Lambda(p,p);
arma::colvec one(p,arma::fill::ones);
arma::mat Cmatrix(p,p);

Cmatrix=Rcpp::as<arma::mat>(C_matrix);



for (int i =0;i <p;i++){
        for (int j =0;j<p;j++){
            betam(i,j) = fabs(beta[i] + C_matrix(i,j)*beta[j]);

        }
    }

count = -1;
for (int i=0;i<p;i++){
    for (int j = 0;j<p;j++){
        if ( i > j)
        {
            count = count + 1;
            b_mu_ig[count]     = (b*sqrt(sigma2))/(fabs(betam(i,j)*sqrt(r))); //betam is already nonnegative
            b_lambda_ig[count]=  pow(b,2);

           }
           if(i==j){
            a_mu_ig[i] = a*sqrt(sigma2)/(sqrt(r)*fabs(beta[i]));
            a_lambda_ig[i] =pow(a,2);
           }

    }
}

b_ig = inversegauss(d,b_mu_ig,b_lambda_ig);
a_ig = inversegauss(p,a_mu_ig,a_lambda_ig);

count  = -1;

for (int i=0;i<p;i++){
    for (int j = 0;j<p;j++){
        if ( i > j)
        {
            count = count + 1;
            Lambda_off(i,j) =b_ig(count);
           }

}
}

Lambda = Rcpp::as<arma::mat>(Lambda_off);
Lambda = Lambda + Lambda.t();
Lambda = Lambda%Cmatrix ;




Lambda.diag() = arma::sum(arma::abs(Lambda),1) + one + Rcpp::as<arma::colvec>(a_ig);

Lambda = r*Lambda;

Lambda_off = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(Lambda));

return  Lambda_off;

}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double samplesigma2inC(Rcpp::NumericMatrix Lambda,Rcpp::NumericMatrix XtX,Rcpp::NumericVector XtYminusYbar,Rcpp::NumericVector Y){
int n = Y.size();
int p = Lambda.ncol();
double shape = n/2;
arma::vec rate(1);

arma::vec y(n) ;
y = Rcpp::as<arma::vec>(Y);
arma::colvec xtyminusybar(p);
xtyminusybar = Rcpp::as<arma::colvec>(XtYminusYbar);
arma::mat xtx (p,p);
xtx = Rcpp::as<arma::mat>(XtX);
arma::mat lambda(p,p);
lambda = Rcpp::as<arma::mat>(Lambda);


lambda = xtx+lambda;
lambda = arma::inv(lambda);

rate = arma::sum(arma::square(y))/2 - (xtyminusybar.t())*lambda*xtyminusybar/2;

double rate1;
rate1 = Rcpp::as<double>(Rcpp::wrap(rate));

return 1/(R::rgamma(shape,1/rate1));

}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double samplehyperainC(Rcpp::NumericVector Beta,double sigma2,double r,double ahyper){

//variable declaration
int p =Beta.size() ;
double a_shape = 1.0;
double a_rate =0;// =  r*arma::accu(arma::abs(beta))/sqrt(sigma2)/2;

for(int i=0;i<p;i++){
    a_rate = a_rate + fabs(Beta[i]);
}

a_rate = a_rate*r/sqrt(sigma2)/2 +ahyper;

return R::rgamma(a_shape,1/a_rate); //Rgamma takes shape and scale//
}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double samplehyperbinC(Rcpp::NumericVector Beta,Rcpp::NumericMatrix CMatrix,double sigma2,double r,double bhyper){

int p =Beta.size();
arma::colvec beta(p);
beta     = Rcpp::as<arma::vec>(Beta);
arma::mat C_matrix(p,p);
C_matrix = Rcpp::as<arma::mat>(CMatrix);

arma::colvec col_one(p);
col_one.fill(1.0);


C_matrix.diag() = (-1)*col_one;

double b_shape = 1.0;
arma::mat t(p,p);
arma::mat prod(p,p);

prod =beta*(col_one.t());

// t = arma::abs(prod+C_matrix%prod.t());
for(int i=0;i<p;i++){
    for(int j=0;j<p;j++){
        t(i,j) = fabs(prod(i,j) + C_matrix(i,j)*prod(j,i));
    }
}

double z    = arma::accu(t);
double b_rate  = r*z/4/sqrt(sigma2) + bhyper;


return R::rgamma(b_shape,1/b_rate); //Rgamma takes shape and scale//
}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double samplehyperrinC(Rcpp::NumericVector Beta,double sigma2,Rcpp::NumericMatrix C_Matrix,double a,double b, double rhyper){

double b_r = 0.01;
int p = Beta.size();


double r_shape = p/2.0 + rhyper;


arma::colvec beta(p);
beta = Rcpp::as<arma::colvec>(Beta);

arma::mat C_matrix(p,p);
C_matrix = Rcpp::as<arma::mat>(C_Matrix);


arma::colvec col_one(p,arma::fill::ones);

C_matrix.diag() = (-1)*col_one;

arma::mat prod(p,p);
arma::mat t(p,p);

prod =beta*(col_one.t());

for(int i=0;i<p;i++){
    for(int j=0;j<p;j++){
        t(i,j) = fabs(prod(i,j) + C_matrix(i,j)*prod(j,i));
    }
}

arma::colvec beta_mult = arma::square(beta);

double r_rate = arma::accu(beta_mult)/sigma2/2 + a*arma::accu(arma::abs(beta))/sqrt(sigma2)/2 + b*(arma::accu(t))/4/sqrt(sigma2);

r_rate = r_rate+b_r/2;

//double r_rate = b*(arma::accu(t))/(4*sqrt(sigma2));

double r = R::rgamma(r_shape,1/r_rate);// rgamma takes shape and scale//

return r;
}



//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List BVSG(Rcpp::NumericMatrix X,Rcpp::NumericVector Y,Rcpp::NumericVector Y0, double a_init, double b_init, double r_init, Rcpp::NumericVector beta_init, Rcpp::NumericMatrix Lambda_init, double sigma2_init, double ahyper, double bhyper, double rhyper, int Nburn, int Niter, int Nthin ){

//std::cout << "The error line is after this " << std::endl;

double a = a_init;
double b = b_init;
double r = r_init;

int n = Y.size();
int p = X.ncol();

int N = (Niter)/Nthin;
int ind;

Rcpp::NumericVector ahist(N);
Rcpp::NumericVector bhist(N);
Rcpp::NumericVector rhist(N);
Rcpp::NumericVector sigma2hist(N);
Rcpp::NumericVector sigma2m(N)  ;


arma::mat Betahist(N,p);

Rcpp::NumericVector beta(p);
beta     = beta_init;
Rcpp::NumericMatrix Lambda(p,p);
Lambda   = Lambda_init;
Rcpp::NumericMatrix C_matrix(p,p);
double sigma2   = sigma2_init;

arma::mat x(n,p);
x            = Rcpp::as<arma::mat>(X);
arma::colvec y(n) ;
y            = Rcpp::as<arma::colvec>(Y);

arma::mat xtx(p,p);
xtx          = (x.t())*x;
arma::colvec xtzminuszbar(p);

Rcpp::NumericMatrix XtX(p,p);
XtX           = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(xtx));

Rcpp::NumericVector XtZminusZbar(p);
Rcpp::NumericVector Z(p) ;

arma::mat covmat(p,p) ;


arma::vec Betahat(p);
arma::mat Lambdahat(p,p);

Z = Y;

double mean_sigma2 ;

    for(int i=0;i<Nburn;i++){


        xtzminuszbar = (x.t())*(Rcpp::as<arma::colvec>(Z));
        XtZminusZbar = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(xtzminuszbar));

        C_matrix   = sampleC_matrixinC(beta,sigma2,r,a,b);
        sigma2     = samplesigma2inC(Lambda,XtX,XtZminusZbar,Z);
        beta       = sampleBetainC(sigma2,Lambda,XtX,XtZminusZbar);
        Lambda     = sampleLambdainC(beta,C_matrix,sigma2,r,a,b);

        a         = samplehyperainC(beta,sigma2,r,ahyper) ;
        b         = samplehyperbinC(beta,C_matrix,sigma2,r,bhyper) ;
        r         = samplehyperrinC(beta,sigma2,C_matrix,a,b,rhyper);

        sigma2m[ i ] = sigma2 ;
        double L  = i - 100 + 1;
        double U  = i ;


        if ( i < 100 )
        {

        Z          = sampleTruncY(X,Y0,beta,sigma2) ;
        } else

        {

        mean_sigma2 = 0;

        for ( int j = L; j <= U ;++ j) {
            mean_sigma2 = mean_sigma2 + sigma2m[j] ;
        }

        mean_sigma2 = mean_sigma2/100;

        Z          = sampleTruncY(X,Y0,beta,mean_sigma2) ;
        }




        if( (1+i) % Nthin == 0)
            {

                ind       =     i/Nthin;
                ahist[ind]      = a;
                bhist[ind]      = b;
                rhist[ind]      = r;
                sigma2hist[ind] = sigma2;

                Betahist.row(ind)  = Rcpp::as<arma::rowvec>(Rcpp::wrap(beta));
            }

    }

    for(int i=0;i<Niter-Nburn;i++){

        xtzminuszbar = (x.t())*(Rcpp::as<arma::colvec>(Z));
        XtZminusZbar = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(xtzminuszbar));

        C_matrix   = sampleC_matrixinC(beta,sigma2,r,a,b);
        sigma2     = samplesigma2inC(Lambda,XtX,XtZminusZbar,Z);
        beta       = sampleBetainC(sigma2,Lambda,XtX,XtZminusZbar);

        Lambda    = sampleLambdainC(beta,C_matrix,sigma2,r,a,b);

        a         = samplehyperainC(beta,sigma2,r,ahyper) ;
        b         = samplehyperbinC(beta,C_matrix,sigma2,r,bhyper) ;
        r         = samplehyperrinC(beta,sigma2,C_matrix,a,b,rhyper);

        sigma2m[ i + Nburn ] = sigma2 ;
        double L  = i + Nburn - 100 + 1;
        double U  = i + Nburn ;


        mean_sigma2 = 0;

        for ( int j = L; j <= U ;++ j) {
            mean_sigma2 = mean_sigma2 + sigma2m[j] ;
        }

        mean_sigma2 = mean_sigma2/100;
        Z          = sampleTruncY(X,Y0,beta,mean_sigma2) ;




        if ((Nburn+i)%Nthin == 0)

            {

                ind         = (Nburn+i)/Nthin;

                ahist[ind]      = a;
                bhist[ind]      = b;
                rhist[ind]      = r;
                sigma2hist[ind]   = sigma2;

                Betahist.row(ind)  = Rcpp::as<arma::rowvec>(Rcpp::wrap(beta));
                sigma2hist[ind]   = sigma2;

                Lambdahat             = Lambdahat + Rcpp::as<arma::mat>(Rcpp::wrap(Lambda));
                Betahat               = Betahat   + Rcpp::as<arma::vec>(Rcpp::wrap(beta));

            }


    }


    int total_iter = (Niter-Nburn)/Nthin;

    Lambdahat = Lambdahat/total_iter;
    Betahat   = Betahat/total_iter;


    return Rcpp::List::create(Rcpp::Named("Betahat")=Betahat,Rcpp::Named("Lambdahat")=Lambdahat,Rcpp::Named("ahist")=ahist,Rcpp::Named("bhist")=bhist,Rcpp::Named("rhist")=rhist,Rcpp::Named("sigma2hist")=sigma2hist,Rcpp::Named("Betahist")=Betahist);

}

