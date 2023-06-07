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
           

  int<lower=1> fixedyear; // centering value for years
 
 // spatial neighbourhood information
  int<lower=1> N_edges;
  array [N_edges] int<lower=1, upper=nroutes> node1;  // node1[i] adjacent to node2[i]
  array [N_edges] int<lower=1, upper=nroutes> node2;  // and node1[i] < node2[i]

  // data for spline s(cov) on trend
  int<lower=1> nknots_trend;  // number of knots in the basis function for year
  matrix[nroutes, nknots_trend] trend_pred_basis; // basis function matrix

  // data for spline s(cov2) on alpha
  int<lower=1> nknots_alpha;  // number of knots in the basis function for year
  matrix[nroutes, nknots_alpha] alpha_pred_basis; // basis function matrix


}

parameters {

  vector[nroutes] beta_raw;
  real BETA; 
  vector[nknots_trend] BETA_raw_hab; 

  vector[nroutes] alpha_raw;
  real ALPHA; 
  vector[nknots_alpha] ALPHA_raw_hab; 

  real eta; //first-year intercept
  
  vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // inverse of sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
 real<lower=0> sdbeta;    // sd of slopes 
 real<lower=0> sdbeta_hab;    // sd of slopes 
  real<lower=0> sdalpha;    // sd of intercepts
  real<lower=0> sdalpha_hab;    // sd of intercepts

  
}

transformed parameters{
  
   vector[nroutes] beta;
   vector[nroutes] beta_trend;
   vector[nroutes] SMOOTH_pred_trend; 
   
  vector[nroutes] alpha;
  vector[nroutes] SMOOTH_pred_alpha;
  vector[nroutes] alpha_trend;
  vector[nknots_trend] BETA_hab; 
  vector[nknots_alpha] ALPHA_hab; 

  vector[nobservers] obs;
  real phi;
  vector[ncounts] E;           // log_scale additive likelihood


// covariate effect on intercepts and slopes
 
 
   BETA_hab = sdbeta_hab*BETA_raw_hab;
   SMOOTH_pred_trend = trend_pred_basis * BETA_hab; 
   beta_trend = (sdbeta*beta_raw) + BETA;
  
   ALPHA_hab = sdalpha_hab*ALPHA_raw_hab;
   SMOOTH_pred_alpha = alpha_pred_basis * ALPHA_hab; 
   alpha_trend = (sdalpha*alpha_raw) + ALPHA;
   
   for(s in 1:nroutes){
     alpha[s] = alpha_trend[s] + (SMOOTH_pred_alpha[s]);
     beta[s] = beta_trend[s]  + (SMOOTH_pred_trend[s]);
   }
   

   obs = sdobs*obs_raw;

  phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

  for(i in 1:ncounts){
    E[i] =  alpha[route[i]] + beta[route[i]] * (year[i]-fixedyear)  + obs[observer[i]] + eta*firstyr[i];
  }
  
  
}

model {

    // spatially varying trend
  beta_raw ~ icar_normal(nroutes, node1, node2); //normal(0,1);//random slope effects
  alpha_raw ~ icar_normal(nroutes, node1, node2);

  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance

  sdobs ~ normal(0,0.3); //prior on sd of gam hyperparameters
 
  obs_raw ~ normal(0,1);//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);

  count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  BETA_raw_hab ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA ~ normal(0,1);// prior on fixed effect mean intercept
  ALPHA_raw_hab ~ normal(0,1);
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and slopes by strata
  sdalpha ~ normal(0,1); //prior on sd of intercept variation
  sdalpha_hab ~ normal(0,1); //prior on sd of habitat smooth parameters

  sdbeta ~ normal(0,0.1);//~ normal(0,0.05); prior on sd of trend variation
  sdbeta_hab ~ normal(0,1);//~ normal(0,0.05); // prior on sd of habitat trend smooth



}

 generated quantities {


   array[nroutes,nyears] real<lower=0> nsmooth;
   array[nroutes,nyears] real<lower=0> nsmooth_no_habitat;
    array[nyears] real<lower=0> NSmooth;
    array[nyears] real<lower=0> NSmooth_no_habitat;
    real T; // full end-point trend of slope-based trajectory
    real T_no_habitat; // end-point trend of slope-based trajectory without habitat-component of slope
    real CH; // total change transformation of trend 
    real CH_no_habitat; //total change transformation of trend without habitat-component of slope
    real CH_dif; // difference between CH and CH_no_habitat
    real T_dif; // difference between T and T_no
  real retrans_obs = 0;// 0.5*(sdobs^2);

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



 

