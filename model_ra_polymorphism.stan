/*
  Stan model for modelling polymorphisms of kinetic parameters
  This has been modified from file _m3.stan, by making gamma_drug hierarchical
*/

data {
  //counts
  int  n;
  int  n_drugs;
  int  n_studies;
  int  n_polymorphisms;  
  //identifiers
  int  drug[n];
  int  study[n];
  int  polymorphism[n]; //code the polymorphism with EF=0 as 1,
                        //and  the one with EF=1 as n_polymorphisms
  real cl_vs_auc; //set 1 for CL and -1 for AUC
  real<lower=0> fm_mean_prior[n_drugs];
  real<lower=0> fm_sd_prior[n_drugs];
  
  //observed data on means & variation:
  real ss[n];
  real lgm[n];
  real<lower=0> lv[n];
}

parameters {
  real mu_drug[n_drugs];     
  real mu_study[n_studies];  
  real<lower=0> sigma_study;
  real gamma_group[n_polymorphisms-1];
  real gamma_drug[n_drugs];
  real mean_gamma;
  real<lower=0> sd_gamma;
  
  // random effect of group:
  real<lower=0, upper=1> mean_ef[n_polymorphisms-2];
  real<lower=0, upper=1> fm[n_drugs];
  real<lower=0, upper =1> ef[n];
  real<lower=0> sigma_ef;
}

transformed parameters {
  real logratio[n];
  real mu[n];
  real sigma[n]; 
  for(i in 1:n){
    //for AUC we need -1 before log() because R_AUC = 1/R_CL
    //for CL it's +1
    logratio[i] = log(ef[i]*fm[drug[i]] + 1-fm[drug[i]]);
    mu[i] = mu_drug[drug[i]] + mu_study[study[i]] + cl_vs_auc*logratio[i];
    sigma[i] = exp(gamma_drug[drug[i]]);
    if(polymorphism[i] > 1)
      sigma[i] = sigma[i]*exp(gamma_group[polymorphism[i]-1]);
  }
}

model {
  //priors dealing with means (mu):
  mu_drug  ~ normal(0, 10);
  mu_study ~ normal(0, sigma_study);
  sigma_study ~ normal(0, 5);
  
  //priors for fm:
  for(i in 1:n_drugs)
    fm[i] ~ normal(fm_mean_prior[i], fm_sd_prior[i]);
  sigma_ef ~ normal(0, 1);
  
  //priors dealing with variances (gamma):
  mean_gamma ~ normal(0, 2.5);
  sd_gamma ~ normal(0, 2.5);
  gamma_drug ~ normal(mean_gamma, sd_gamma);
  gamma_group ~ normal(0, 5);
  for(i in 1:(n_polymorphisms-2))
    mean_ef[i] ~ uniform(0,1);
  
  for(i in 1:n) {
    //enzyme function:
    if(polymorphism[i] == 1) 
      ef[i] ~ normal(0, sigma_ef);
    if(polymorphism[i] == n_polymorphisms) 
      ef[i] ~ normal(1, sigma_ef);
    if((polymorphism[i] > 1) && (polymorphism[i] < n_polymorphisms))
      ef[i] ~ normal(mean_ef[polymorphism[i]-1], sigma_ef);
    
    //observed quantities:
    lgm[i] ~ normal(mu[i], sigma[i]/sqrt(ss[i]));
    //chi-sq is Gamma(njk - 1 / 2, tau_j*(n_jk-1)/2)
    lv[i]  ~ gamma((ss[i] - 1) / 2, (ss[i]-1)/(2*(sigma[i]^2))); 
  }
}
