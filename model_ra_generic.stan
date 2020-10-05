data {
  int  n;
  int  n_drugs;
  int  n_studies;
  int  n_groups;  
  int  drug[n];
  int  study[n];
  int  group[n];
  
  //observed data on means & variation:
  real ss[n];
  real lgm[n];
  real<lower=0> lv[n];
}

parameters {
  real mu_drug[n_drugs];     
  real mu_study[n_studies];  
  real<lower=0> sigma_study;
  real gamma_group[n_groups];
  real gamma_drug[n_drugs];
  
  // random effect of group:
  // real logratio[n];
  // real mu_group[n_groups - 1];
  // real<lower=0> sigma_group;
  
  // fixed effect of group:
  real logratio[n_groups - 1];
}

model {
  // Both true mu and sigma are simply introduced to clean 
  // up notation, rather than as extra parameters to output.
  // Hence they are defined here, in the model block.
  real sigma[n]; 
  real mu[n];
  for(i in 1:n){
    if(group[i] != 1)
      // random effect
      // mu[i] = mu_drug[drug[i]] + mu_study[study[i]] + logratio[n];
      // fixed effect
      mu[i] = mu_drug[drug[i]] + mu_study[study[i]] + logratio[group[i] - 1];
    else
      mu[i] = mu_drug[drug[i]] + mu_study[study[i]];
    sigma[i] = exp(gamma_drug[drug[i]] + gamma_group[group[i]]);
  }
  
  //priors dealing with means (mu):
  mu_drug ~ normal(0, 100);
  for(i in 1:n_studies)
    mu_study[i] ~ normal(0, sigma_study);
  sigma_study ~ normal(0, 5);
  gamma_drug ~ normal(0, 10);
  gamma_group ~ normal(0, 10);
  
  // random effect of group:
  // mu_group ~ normal(0, 100);
  // sigma_group ~ normal(0, 5);
  
  // fixed effect of group:
  logratio ~ normal(0, 5);
  
  for(i in 1:n){
    // random effect of group:
    // if(group[i] == 1)
    //   logratio[i] ~ normal(0, .0001);
    // else
    //   logratio[i] ~ normal(mu_group[group[i] - 1], sigma_group);
    //observed quantities:
    lgm[i] ~ normal(mu[i], sigma[i]/sqrt(ss[i]));
    //chi-sq is Gamma(njk - 1 / 2, tau_j*(n_jk-1)/2)
    lv[i]  ~ gamma((ss[i] - 1) / 2, (ss[i]-1)/(2*(sigma[i]^2))); 
  }
}
