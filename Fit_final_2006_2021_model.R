## Fitting the habitat-change model with iCAR trend component to 2006 - 2021 BBS data
setwd("C:/GitHub/Jeffreys_etal") #useful if running in standaolone R session (avoids Rstudio fails)
setwd("C:/Users/SmithAC/Documents/GitHub/Jefferys_etal")

library(tidyverse)
library(cmdstanr)
library(patchwork)
library(sf)




output_dir <- "output"
species <- "Rufous Hummingbird" 


species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)

spp1 <- "habitat"

spp <- paste0("_",spp1,"_")


for(firstYear in c(1985,2006)){
  
  lastYear <- ifelse(firstYear == 1985,2005,2021)
  




out_base <- paste0(species_f,spp,firstYear,"_",lastYear)




sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")


load(sp_data_file)

# trend habitat effects are not changed, but the intercept effect is
# removes the optional spatial components for intercepts in 1985-2005
stan_data[["fit_spatial"]] <- ifelse(firstYear == 1985,0,1)

   mod.file = paste0("models/slope",spp,"route_NB.stan")


slope_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))

stanfit <- slope_model$sample(
  data=stan_data,
  refresh=400,
  iter_sampling=2000,
  iter_warmup=2000,
  parallel_chains = 4)

summ <- stanfit$summary()
print(paste(species, stanfit$time()[["total"]]))

saveRDS(stanfit,
        paste0(output_dir,"/",out_base,"_stanfit.rds"))

saveRDS(summ,
        paste0(output_dir,"/",out_base,"_summ_fit.rds"))

summ %>% arrange(-rhat)
summ %>% filter(variable %in% c("BETA","BETA_hab"))
summ %>% filter(grepl("T",variable))
summ %>% filter(grepl("CH",variable))



}


# graphing ----------------------------------------------------------------


firstYear <- 2006
lastYear <- 2021

out_base <- paste0(species_f,spp,firstYear,"_",lastYear)

sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")

load(sp_data_file)

summ <- readRDS(paste0(output_dir,"/",out_base,"_summ_fit.rds"))


mn0 <- new_data %>% 
  group_by(routeF) %>% 
  summarise(mn = mean(count),
            mx = max(count),
            ny = n(),
            fy = min(year),
            ly = max(year),
            sp = max(year)-min(year))

route_map_2006 <- route_map 

exp_t <- function(x){
  y <- (exp(x)-1)*100
}


# plot trends -------------------------------------------------------------


base_strata_map <- bbsBayes2::load_map("bbs_usgs")


strata_bounds <- st_union(route_map) #union to provide a simple border of the realised strata
bb = st_bbox(strata_bounds)
xlms = as.numeric(c(bb$xmin,bb$xmax))
ylms = as.numeric(c(bb$ymin,bb$ymax))

betas1 <- summ %>% 
  filter(grepl("beta[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp_t(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "full") %>% 
  select(routeF,mean,sd,parameter) %>% 
  rename(trend = mean,
         trend_se = sd)

alpha1 <- summ %>% 
  filter(grepl("alpha[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}"))) %>% 
  select(routeF,median,sd) %>% 
  rename(abundance = median,
         abundance_se = sd)
betas1 <- betas1 %>% 
  inner_join(.,alpha1)

betas2 <- summ %>% 
  filter(grepl("beta_trend[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp_t(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "no habitat") %>% 
  select(routeF,mean,sd,parameter) %>% 
  rename(trend = mean,
         trend_se = sd)
betas2 <- betas2 %>% 
  inner_join(.,alpha1,by = "routeF") %>% 
  inner_join(.,mn0,by = "routeF")


betas <- bind_rows(betas1,betas2)

plot_map <- route_map_2006 %>% 
  left_join(.,betas,
            by = "routeF",
            multiple = "all") 

breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
lgnd_head <- "Mean Trend\n"
trend_title <- "Mean Trend"
labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
labls = paste0(labls, " %/year")
plot_map$Tplot <- cut(plot_map$trend,breaks = c(-Inf, breaks, Inf),labels = labls)


map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                 "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
names(map_palette) <- labls



map <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = plot_map,
          aes(colour = Tplot,
              size = abundance))+
  scale_size_continuous(range = c(0.05,2),
                        name = "Mean Count")+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0(lgnd_head))+
  coord_sf(xlim = xlms,ylim = ylms)+
  theme_bw()+
  labs(title = paste(firstYear,"-",lastYear))+
  facet_wrap(vars(parameter))



#map

# pdf(paste0("Figures/Four_trends_model_comparison_",species_f,".pdf"),
#     height = 8,
#     width = 8)
# print(map)
# dev.off()

map_se <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = plot_map,
          aes(colour = trend_se,
              size = abundance_se))+
  scale_size_continuous(range = c(0.05,2),
                        name = "SE of Mean Count",
                        trans = "reverse")+
  scale_colour_viridis_c(aesthetics = c("colour"),
                         guide = guide_legend(reverse=TRUE),
                         name = paste0("SE of Trend"))+
  coord_sf(xlim = xlms,ylim = ylms)+
  theme_bw()+
  labs(title = paste(firstYear,"-",lastYear))+
  facet_wrap(vars(parameter))




#print(map2 / map_se2)

pdf(paste0("Figures/Figure_supplement_1_Trend_map_w_habitat_and_withing_",species_f,".pdf"),
    height = 10.5,
    width = 7.5)


print(map / map_se + plot_layout(guides = "collect"))


dev.off()


pdf(paste0("Figures/Figure_4trend.pdf"),
    height = 5,
    width = 7)


print(map)


dev.off()

  
  
  
  