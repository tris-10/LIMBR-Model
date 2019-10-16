data{
  // get the data set
  int<lower=0> N;   // number of domain/exon level
  int<lower=0> J;   // number of gene level
  int <lower=1,upper=J> gene[N];  //index of gene, restricted from 1 to J
  int <lower=1,upper=N> domain[N];  //index of domains, restricted from 1 to N
  vector[N] xij;   //x at domain level
  vector[N] yij; //y at domain level
  vector[N] hasLowSubReg_ij; //indicator for low number of sub-regions 
}
parameters{
  // specify the parameter we want to know
  real beta;  //common slope for the domain level
  real mu;      //common intercept for the domain level
  
  vector[N] offset_alpha_ij;
  real <lower=0> sigma_aj_Sqrd[J];  //variance of intercept at domain level
  vector[J] alpha_j; //random intercept for the gene level
  real <lower=0.00001> sigma_a_Sqrd;  //variance of intercept at gene level
  real <lower=0> sigma_Sqrd; //variance of yij
  real <lower=0> eps_aj; //hyperparameter for sigma_aj_Sqrd
  real <lower=0> eps_a; //hyperparameter for sigma_a
  real <lower=0> eps; //hyperparameter for sigma
}
transformed parameters{
  // combine the alpha terms 
  vector[N] combinedAlphas;
  vector[N] alpha_ij; //random intercept for the domain level
  for (i in 1:N){
    if(hasLowSubReg_ij[i] == 1){
      alpha_ij[i] = 0; //offset_alpha_ij[i] * sqrt(tau);
    }
    else{
      alpha_ij[i] = offset_alpha_ij[i] * sqrt( sigma_aj_Sqrd[gene[i]] ); //tau*
    }
    combinedAlphas[i] = alpha_ij[i]+alpha_j[gene[i]];  
  }
  
}
model {
  // give the prior distribution
  vector[N] mu_ij; #domain level model
  for (i in 1:N)
    mu_ij[i] = mu+beta*xij[i]+alpha_ij[domain[i]]+alpha_j[gene[i]]; //specify the gene group
  beta ~normal(0,50);
  mu~normal(0,50);
  eps ~ uniform(0.00001,100);
  eps_a ~ uniform(0.00001,100);
  eps_aj ~ uniform(0.00001,100);
  sigma_Sqrd ~ inv_gamma(eps,eps);//uniform(0.00001, 20);
  sigma_a_Sqrd ~ inv_gamma(eps_a,eps_a);//uniform(0.00001,10);
  sigma_aj_Sqrd ~ inv_gamma(eps_aj,eps_aj);
  alpha_j ~ normal(0, sqrt(sigma_a_Sqrd));
  for (i in 1:N){
    offset_alpha_ij[i] ~ normal(0,1);
  }
  yij ~ normal(mu_ij,sqrt(sigma_Sqrd));
}
