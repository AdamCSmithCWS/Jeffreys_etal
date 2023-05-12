

library(tidyverse)




output_dir <- "output"



firstYear <- 1985
lastYear <- 2021


species <- "Rufous Hummingbird" 


species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)


sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")


load(sp_data_file)


summ <- NULL
  
  
for(spp1 in c("habitat","nonhabitat")){
  
  
  spp <- paste0("_",spp1,"_")
  
  
  out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
  


tmp <- readRDS(paste0(output_dir,"/",out_base,"_summ_fit.rds")) %>% 
  mutate(model = spp1)

summ <- bind_rows(summ,
                  tmp)



}


BETA <- summ %>% 
  filter(variable == "BETA")


betas <- summ %>% 
  filter(grepl("beta[",
               variable,
               fixed = TRUE))

alphas <- summ %>% 
  filter(grepl("alpha[",
               variable,
               fixed = TRUE))


alpha_habs <- summ %>% 
  filter(grepl("alpha_hab[",
               variable,
               fixed = TRUE))

BETA_habs <- summ %>% 
  filter(variable == "BETA_hab")


beta_habs <- summ %>% 
  filter(grepl("beta_hab[",
               variable,
               fixed = TRUE))

ALPHA_habs <- summ %>% 
  filter(variable == "ALPHA_hab")




