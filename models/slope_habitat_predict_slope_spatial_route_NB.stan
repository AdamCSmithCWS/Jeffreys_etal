// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends
// and no random year-effects - slope only

//iCAR function
 functions {
   real icar_normal_lpdf(vector bb, int ns, array[] int n1, array[] int n2) {
     return -0.5 * dot_self(bb[n1] - bb[n2])
       + normal_lpdf(sum(bb) | 0, 0.001 * ns); //soft sum to zero constraint on bb
  }
 }


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

  array [nroutes] real route_habitat_slope;              

  int<lower=1> fixedyear; // centering value for years
 
 // spatial neighbourhood information
  int<lower=1> N_edges;
  array [N_edges] int<lower=1, upper=nroutes> node1;  // node1[i] adjacent to node2[i]
  array [N_edges] int<lower=1, upper=nroutes> node2;  // and node1[i] < node2[i]


}

parameters {

  vector[nroutes] beta_raw;
  real BETA; 
  vector[nroutes] beta_raw_hab;
  real BETA_hab; 

  vector[nroutes] alpha_raw;
  //vector[nroutes] alpha_raw_hab;
  real ALPHA; 
  real ALPHA_hab; 

  real eta; //first-year intercept
  
  vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // inverse of sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
 real<lower=0> sdbeta;    // sd of slopes 
 real<lower=0> sdbeta_hab;    // sd of slopes 
  real<lower=0> sdalpha;    // sd of intercepts
  //real<lower=0> sdalpha_hab;    // sd of intercepts

  
}

transformed parameters{
  
   vector[nroutes] beta;
   vector[nroutes] beta_trend;
   vector[nroutes] beta_hab;
   //vector[nroutes] beta_dif;
   //real BETA_dif;
   
  vector[nroutes] alpha;
  //vector[nroutes] alpha_hab;
  vector[nobservers] obs;
  real phi;
  vector[ncounts] E;           // log_scale additive likelihood


// covariate effect on intercepts and slopes
 
   beta_trend = (sdbeta*beta_raw) + BETA;
   beta_hab = (sdbeta_hab*beta_raw_hab) + BETA_hab;
   
   //alpha_hab = (sdalpha_hab*alpha_raw_hab) + ALPHA_hab;
   
   for(s in 1:nroutes){
     alpha[s] = (sdalpha*alpha_raw[s]) + ALPHA  + (ALPHA_hab*route_habitat[s]);
     beta[s] = beta_trend[s]  + (beta_hab[s]*route_habitat_slope[s]);
     //beta_dif[s] = beta[s]-beta_trend[s];
   }
   
  // BETA_dif = mean(beta_dif);
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
  
  // spatially varying coefficient of habitat change on trend
  beta_raw_hab ~ icar_normal(nroutes, node1, node2);//~ normal(0,1);//habitat slope effects
  //sum(beta_raw_hab) ~ normal(0,0.001*nroutes);
  alpha_raw ~ normal(0,1);
  sum(alpha_raw) ~ normal(0,0.001*nroutes);
  // alpha_raw_hab ~ normal(0,1);
  // sum(alpha_raw_hab) ~ normal(0,0.001*nroutes);

  
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
  //sdalpha_hab ~ normal(0,2); //prior on sd of intercept variation

  sdbeta ~ normal(0,0.1);//~ normal(0,0.05); //boundary avoiding prior on sd of slope spatial variation w mean = 0.04 and 99% < 0.13
  sdbeta_hab ~ normal(0,0.1);//~ normal(0,0.05); //boundary avoiding prior on sd of slope spatial variation w mean = 0.04 and 99% < 0.13



}

 generated quantities {


   array[nroutes,nyears] real<lower=0> nsmooth;
   array[nroutes,nyears] real<lower=0> nsmooth_no_habitat;
    array[nyears] real<lower=0> NSmooth;
    array[nyears] real<lower=0> NSmooth_no_habitat;
    real T;
    real T_no_habitat;
    real CH;
    real CH_no_habitat;
    real CH_dif;
    real T_dif;
  real retrans_obs = 0.5*(sdobs^2);

// intercepts and slopes


 for(y in 1:nyears){

      for(s in 1:nroutes){



      nsmooth[s,y] = exp(alpha[s] + beta[s] * (y-fixedyear)  + retrans_obs);//
      nsmooth_no_habitat[s,y] = exp(alpha[s] + beta_trend[s] * (y-fixedyear)  + retrans_obs);//
        }
  NSmooth[y] = mean(nsmooth[,y]);
  NSmooth_no_habitat[y] = mean(nsmooth_no_habitat[,y]);
  
    }
    T = 100*(((NSmooth[nyears]/NSmooth[1])^(1.0/nyears))-1);
    
    
    T_no_habitat = 100*(((NSmooth_no_habitat[nyears]/NSmooth_no_habitat[1])^(1.0/nyears))-1);
    
        CH = 100*((NSmooth[nyears]/NSmooth[1])-1);
    
    
    CH_no_habitat = 100*((NSmooth_no_habitat[nyears]/NSmooth_no_habitat[1])-1);
    
    CH_dif = CH-CH_no_habitat;
    T_dif = T-T_no_habitat;
    
  }



 
