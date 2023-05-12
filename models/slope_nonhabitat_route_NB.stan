// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends
// and no random year-effects - slope only



data {
  int<lower=1> nroutes;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  int<lower=1> nobservers;
 
  array [ncounts] int<lower=0> count;              // count observations
  array [ncounts] int<lower=1> year; // year index
  array [ncounts] int<lower=1> route; // route index
  array [ncounts] int<lower=0> firstyr; // first year index
  array [ncounts] int<lower=1> observer;              // observer indicators
  array [nroutes] real route_habitat;              

  int<lower=1> fixedyear; // centering value for years
 


}

parameters {

  vector[nroutes] beta_raw;
  real BETA; 

  vector[nroutes] alpha_raw;
  vector[nroutes] alpha_raw_hab;
  real ALPHA; 
  real ALPHA_hab; 

  real eta; //first-year intercept
  
  vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // inverse of sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
 real<lower=0> sdbeta;    // sd of slopes 
  real<lower=0> sdalpha;    // sd of intercepts
  real<lower=0> sdalpha_hab;    // sd of intercepts

  
}

transformed parameters{
  
   vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nroutes] alpha_hab;
  vector[nobservers] obs;
  real phi;
  vector[ncounts] E;           // log_scale additive likelihood


// covariate effect on intercepts and slopes
 
   beta = (sdbeta*beta_raw) + BETA;
   alpha_hab = (sdalpha_hab*alpha_raw_hab) + ALPHA_hab;
   
   for(s in 1:nroutes){
     alpha[s] = (sdalpha*alpha_raw[s]) + ALPHA  + (alpha_hab[s]*route_habitat[s]);
   }
 //  noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;

  phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

  for(i in 1:ncounts){
    E[i] =  alpha[route[i]] + beta[route[i]] * (year[i]-fixedyear)  + obs[observer[i]] + eta*firstyr[i];
  }
  
  
}

model {



  
  beta_raw ~ normal(0,1);//random slope effects
  sum(beta_raw) ~ normal(0,0.001*nroutes);
  alpha_raw ~ normal(0,1);
  sum(alpha_raw) ~ normal(0,0.001*nroutes);
  alpha_raw_hab ~ normal(0,1);
  sum(alpha_raw_hab) ~ normal(0,0.001*nroutes);

  
  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance


  sdobs ~ std_normal(); //prior on sd of gam hyperparameters
 
  obs_raw ~ normal(0,1);//observer effects
 sum(obs_raw) ~ normal(0,0.001*nobservers);

  count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.2);// prior on fixed effect mean slope
  ALPHA ~ student_t(3,0,1);;// prior on fixed effect mean intercept
  ALPHA_hab ~ student_t(10,0,3);;// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and slopes by strata
  sdalpha ~ normal(0,2); //prior on sd of intercept variation
  sdalpha_hab ~ normal(0,2); //prior on sd of intercept variation

  sdbeta ~ normal(0,0.1);//~ normal(0,0.05); //boundary avoiding prior on sd of slope spatial variation w mean = 0.04 and 99% < 0.13



}

 generated quantities {


   array[nroutes,nyears] real<lower=0> nsmooth;
  real retrans_obs = 0.5*(sdobs^2);

// intercepts and slopes
 
   
 for(y in 1:nyears){

      for(s in 1:nroutes){

 

      nsmooth[s,y] = exp(alpha[s] + beta[s] * (y-fixedyear)  + retrans_obs);// 
        }
 
    }
  }



 

