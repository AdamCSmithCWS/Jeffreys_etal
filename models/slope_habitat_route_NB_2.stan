// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends
// and no random year-effects - slope only



data {
  int<lower=1> nroutes;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  int<lower=0> nobservers;
 
  array [ncounts] int<lower=0> count;              // count observations
  array [ncounts] int<lower=1> year; // year index
  array [ncounts] int<lower=1> route; // route index
  array [ncounts] int<lower=0> firstyr; // first year index
  array [ncounts] int<lower=1> observer;              // observer indicators
  array [ncounts] real c_habitat;              
  array [nroutes] real route_habitat;              

   int<lower=1> fixedyear; // centering value for years
  // 

}

parameters {

//  vector[nroutes] beta_raw_hab;
  vector[nroutes] beta_raw;
  
  real BETA; 
 // real BETA_hab; 

  vector[nroutes] alpha_raw;
  vector[nroutes] alpha_raw_hab;
  real ALPHA; 
  real ALPHA_hab; 

  real eta; //first-year intercept
  
  vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // inverse of sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
//  real<lower=0> sdbeta_hab;    // sd of slopes 
  real<lower=0> sdbeta;    // sd of slopes 
// real<lower=0> sdbeta_rand;    // sd of slopes 
  real<lower=0> sdalpha;    // sd of intercepts
  real<lower=0> sdalpha_hab;    // sd of intercepts

  
}


model {


  vector[ncounts] E;           // log_scale additive likelihood
   //vector[nroutes] beta_rand;
//  vector[nroutes] beta_hab;
 vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nroutes] alpha_hab;
  vector[nobservers] obs;
  real phi;

// covariate effect on intercepts and slopes
   beta_hab = (sdbeta_hab*beta_raw_hab) + BETA_hab;
   //beta_rand = (sdbeta_rand*beta_raw_rand);
   
   //beta = (sdbeta*beta_raw) + BETA;
   
   alpha = (sdalpha*alpha_raw) + ALPHA;
   alpha_hab = (sdalpha_hab*alpha_raw_hab) + ALPHA_hab;
 //  noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;

  for(i in 1:ncounts){
    E[i] =  alpha[route[i]] + beta_hab[route[i]] * (c_habitat[i]) + alpha_hab[route[i]]*route_habitat[route[i]] + obs[observer[i]] + eta*firstyr[i];
 //   E[i] =  alpha[route[i]] + beta[route[i]] * (year[i]-fixedyear) + beta_hab[route[i]] * (c_habitat[i]) + alpha_hab[route[i]]*route_habitat[route[i]] + obs[observer[i]] + eta*firstyr[i];
  }
  
  
  // beta_raw_rand ~ normal(0,1);//random slope effects
  // sum(beta_raw_rand) ~ normal(0,0.001*nroutes);

  
  sdnoise ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance

  phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

  sdobs ~ std_normal(); //prior on sd of gam hyperparameters
 
  obs_raw ~ normal(0,1);//observer effects
//  sum(obs_raw) ~ normal(0,0.001*nobservers);

  count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
//  BETA_hab ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA ~ student_t(10,0,3);;// prior on fixed effect mean intercept
  ALPHA_hab ~ student_t(10,0,3);;// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and slopes by strata
  sdalpha ~ normal(0,2); //prior on sd of intercept variation
  sdalpha_hab ~ student_t(4,0,2); //prior on sd of intercept variation
  // sdbeta_hab ~ gamma(2,50);//~ normal(0,0.05); //boundary avoiding prior on sd of slope spatial variation w mean = 0.04 and 99% < 0.13
  // sdbeta_rand  ~ gamma(2,50);//~ normal(0,0.05); //boundary avoiding prior on sd of slope random variation

 // sdbeta_hab ~ student_t(50,0,0.1);//~ normal(0,0.05); //boundary avoiding prior on sd of slope spatial variation w mean = 0.04 and 99% < 0.13
  sdbeta ~ student_t(50,0,0.1);//~ normal(0,0.05); //boundary avoiding prior on sd of slope spatial variation w mean = 0.04 and 99% < 0.13
  //sdbeta_rand  ~ student_t(10,0,1);//~ normal(0,0.05); //boundary avoiding prior on sd of slope random variation
//  beta_raw_hab ~ normal(0,1);
  beta_raw ~ normal(0,1);
//  sum(beta_raw_hab) ~ normal(0,0.001*nroutes);
  sum(beta_raw) ~ normal(0,0.001*nroutes);
  
  alpha_raw ~ normal(0,1);
  sum(alpha_raw) ~ normal(0,0.001*nroutes);
  alpha_raw_hab ~ normal(0,1);
  sum(alpha_raw_hab) ~ normal(0,0.001*nroutes);


}

 generated quantities {

   //vector[nroutes] beta_rand;
  vector[nroutes] beta_hab;
  vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nroutes] alpha_hab;

// intercepts and slopes
   //beta_hab = (sdbeta_hab*beta_raw_hab);
   //beta_rand = (sdbeta_rand*beta_raw_rand);
   
   beta_hab = (sdbeta_hab*beta_raw_hab) + BETA_hab;
   beta = (sdbeta*beta_raw) + BETA;
   
   alpha = (sdalpha*alpha_raw) + ALPHA;
   alpha_hab = (sdalpha_hab*alpha_raw_hab) + ALPHA_hab;



 }

