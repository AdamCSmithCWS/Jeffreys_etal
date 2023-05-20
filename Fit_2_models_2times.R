## Fitting the BYM model to 1995 - 2021 BBS data
## script currently written to fit the model then save the Stan output to a directory
## 
#setwd("C:/GitHub/iCAR_route_2021")
setwd("C:/Users/SmithAC/Documents/GitHub/Jefferys_etal")
library(tidyverse)
library(cmdstanr)




output_dir <- "output"


for(firstYear in c(1985,2006)){
  
  lastYear <- ifelse(firstYear == 1985,2005,2021)
  

species <- "Rufous Hummingbird" 

  
species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)

for(spp1 in c("habitat_predict_slope","habitat_predict_slope_spatial")[2]){#},"iCAR")){


  spp <- paste0("_",spp1,"_")
  

out_base <- paste0(species_f,spp,firstYear,"_",lastYear)




sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")


load(sp_data_file)


if(spp1 == "nonspatial"){
  
  stan_data[["N_edges"]] <- NULL
  stan_data[["node1"]] <- NULL
  stan_data[["node2"]] <- NULL

  stan_data[["route_habitat_slope"]] <- NULL
  
}

if(spp1 == "habitat_predict_slope"){
  
  stan_data[["N_edges"]] <- NULL
  stan_data[["node1"]] <- NULL
  stan_data[["node2"]] <- NULL
  
 
}




   mod.file = paste0("models/slope",spp,"route_NB.stan")


slope_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))

stanfit <- slope_model$sample(
  data=stan_data,
  refresh=200,
  chains=4, iter_sampling=2000,
  iter_warmup=2000,
  parallel_chains = 4,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 10)

summ <- stanfit$summary()
print(paste(species, stanfit$time()[["total"]]))

saveRDS(stanfit,
        paste0(output_dir,"/",out_base,"_stanfit.rds"))

saveRDS(summ,
        paste0(output_dir,"/",out_base,"_summ_fit.rds"))


} #end models loop


}# end year loop
