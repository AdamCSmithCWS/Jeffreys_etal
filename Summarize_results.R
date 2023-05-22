

library(tidyverse)
library(patchwork)
library(sf)




output_dir <- "output"

summ <- NULL

mean_obs <- NULL

for(firstYear in c(2006,1985)){
  
  
  lastYear <- ifelse(firstYear == 1985,2005,2021)
  

species <- "Rufous Hummingbird" 


species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)


sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")


load(sp_data_file)


  
for(spp1 in c("nonspatial","habitat_predict_slope_spatial")){
  
  
  spp <- paste0("_",spp1,"_")
  
  
  out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
  


tmp <- readRDS(paste0(output_dir,"/",out_base,"_summ_fit.rds")) %>% 
  mutate(model = spp1,
         first_year = firstYear)

summ <- bind_rows(summ,
                  tmp)



}

mean_obst <- new_data %>% 
  group_by(r_year) %>% 
  summarize(mean_obs = mean(count),
            n_routes = n()) %>% 
  mutate(first_year = firstYear)

mean_obs <- bind_rows(mean_obs,mean_obst)

if(firstYear == 1985){
route_map_1985 <- route_map 
}else{
route_map_2006 <- route_map 
}





}#end firstYear

exp_t <- function(x){
  y <- (exp(x)-1)*100
}

mod_sel <- "habitat_predict_slope_spatial"
# year_sel <- 1985
# 
# 
# 
# sp_data_file <- paste0("Data/",species_f,"_",year_sel,"_",lastYear,"_stan_data.RData")
# 
# 
# load(sp_data_file)


BETA <- summ %>% 
  filter(variable == "BETA") %>% 
  mutate(across(2:7,~exp_t(.x)))


betas <- summ %>% 
  filter(grepl("beta[",
               variable,
               fixed = TRUE),
         model == mod_sel)%>% 
  mutate(across(2:7,~exp_t(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")))

beta_trends <- summ %>% 
  filter(grepl("beta_trend[",
               variable,
               fixed = TRUE))%>% 
  mutate(across(2:7,~exp_t(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")))


alphas <- summ %>% 
  filter(grepl("alpha[",
               variable,
               fixed = TRUE),
         model == mod_sel) %>% 
  mutate(routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")))



alpha_habs <- summ %>% 
  filter(grepl("alpha_hab[",
               variable,
               fixed = TRUE),
         model == mod_sel) %>% 
  mutate(routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")))

BETA_habs <- summ %>% 
  filter(variable == "BETA_hab") %>% 
  mutate(across(2:7,~exp_t(.x)))


beta_habs <- summ %>% 
  filter(grepl("beta_hab[",
               variable,
               fixed = TRUE),
         model == mod_sel)%>% 
  mutate(routeF = as.integer(str_extract(pattern = "[[:digit:]]{1,}",
                                         variable)),
         across(2:7,~exp_t(.x)))

ALPHA_habs <- summ %>% 
  filter(variable == "ALPHA_hab")

sds <- summ %>% 
  filter(grepl("sd",
               variable))



N <- summ %>% 
  filter(grepl("NSmooth",
               variable),
         model == mod_sel) %>% 
  mutate(year = as.integer(str_extract(pattern = "[[:digit:]]{1,}",
                            variable)),
         type = ifelse(grepl("_no_habitat",variable,
                             fixed = TRUE),
                       "Habitat Adjusted",
                       "Full"),
         Year = ifelse(first_year == 1985,
                       year + 1984,
                       year + 2005))




traj <- ggplot()+
  geom_ribbon(data = N,
              aes(x = Year,
                  y = median,
                  ymin = q5,
                  ymax = q95,
                  fill = type),
              alpha = 0.2)+
  geom_line(data = N,
            aes(x = Year,
                y = median,
                colour = type))+
  scale_y_continuous(limits = c(0,NA))+
  facet_grid(vars(first_year),
            scales = "free_y")#+
  # geom_point(data = mean_obs,
  #            aes(x = r_year,
  #                y = mean_obs,
  #                alpha = n_routes))
traj





beta_habs1 <- beta_habs %>% 
  filter(first_year == 1985)
t_map1 <- route_map_1985 %>% 
  left_join(.,beta_habs1,
            by = "routeF")
beta_habs2 <- beta_habs %>% 
  filter(first_year == 2006)
t_map2 <- route_map_2006 %>% 
  left_join(.,beta_habs2,
            by = "routeF")


  t_plot1 <- ggplot(data = t_map1)+
  geom_sf(aes(colour = mean))+
  colorspace::scale_color_continuous_diverging(mid = 0,rev = TRUE,
                                               guide = guide_colourbar(title = "Effect of habitat \n change on \n population trend \n 1985-2005",
                                                                    reverse = FALSE))
  
  t_plot2 <- ggplot(data = t_map2)+
    geom_sf(aes(colour = mean))+
    colorspace::scale_color_continuous_diverging(mid = 0,rev = TRUE,
                                                 guide = guide_colourbar(title = "Effect of habitat \n change on \n population trend \n 2006-2021",
                                                                         reverse = FALSE))
  

  
  print(t_plot1 + t_plot2)



# plot trends -------------------------------------------------------------


  base_strata_map <- bbsBayes2::load_map("bbs_usgs")
  
  
  strata_bounds <- st_union(route_map_1985) #union to provide a simple border of the realised strata
  bb = st_bbox(strata_bounds)
  xlms = as.numeric(c(bb$xmin,bb$xmax))
  ylms = as.numeric(c(bb$ymin,bb$ymax))
  
  betas1 <- betas %>% 
    filter(first_year == 1985) %>% 
    select(first_year,routeF,mean,sd) %>% 
    rename(trend = mean,
           trend_se = sd)
  alpha1 <- alphas %>% 
    filter(first_year == 1985) %>% 
    select(first_year,routeF,mean,sd) %>% 
    mutate(abundance = exp(mean),
           abundance_se = sd)
  betas1 <- betas1 %>% 
    inner_join(.,alpha1)
  
  
  plot_map <- route_map_1985 %>% 
    left_join(.,betas1,
              by = "routeF") 
  
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
    labs(title = paste(firstYear,"-",lastYear))
  
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
    labs(title = paste(firstYear,"-",lastYear))
  
    
  
  
  
  strata_bounds <- st_union(route_map_2006) #union to provide a simple border of the realised strata
  bb = st_bbox(strata_bounds)
  xlms = as.numeric(c(bb$xmin,bb$xmax))
  ylms = as.numeric(c(bb$ymin,bb$ymax))
  
  betas2 <- betas %>% 
    filter(first_year == 2006) %>% 
    select(first_year,routeF,mean,sd) %>% 
    rename(trend = mean,
           trend_se = sd)
  alpha2 <- alphas %>% 
    filter(first_year == 2006) %>% 
    select(first_year,routeF,mean,sd) %>% 
    mutate(abundance = exp(mean),
           abundance_se = sd)
  betas2 <- betas2 %>% 
    inner_join(.,alpha2)
  
  
  plot_map2 <- route_map_2006 %>% 
    left_join(.,betas2,
              by = "routeF") 
  
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  lgnd_head <- "Mean Trend\n"
  trend_title <- "Mean Trend"
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %/year")
  plot_map2$Tplot <- cut(plot_map2$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
  
  
  map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                   "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
  names(map_palette) <- labls
  
  
  
  map2 <- ggplot()+
    geom_sf(data = base_strata_map,
            fill = NA,
            colour = grey(0.75))+
    geom_sf(data = plot_map2,
            aes(colour = Tplot,
                size = abundance))+
    scale_size_continuous(range = c(0.05,2),
                          name = "Mean Count")+
    scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                        guide = guide_legend(reverse=TRUE),
                        name = paste0(lgnd_head))+
    coord_sf(xlim = xlms,ylim = ylms)+
    theme_bw()+
    theme(legend.position = "none")+
    labs(title = paste(2006,"-",2021))
  
  # pdf(paste0("Figures/Four_trends_model_comparison_",species_f,".pdf"),
  #     height = 8,
  #     width = 8)
  # print(map)
  # dev.off()
  
  map_se2 <- ggplot()+
    geom_sf(data = base_strata_map,
            fill = NA,
            colour = grey(0.75))+
    geom_sf(data = plot_map2,
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
    theme(legend.position = "none")+
    labs(title = paste(2006,"-",2021))
  
  
  #print(map2 / map_se2)
  
  pdf(paste0("Figures/Trend_map_two_time_periods_",species_f,".pdf"),
  height = 10.5,
  width = 7.5)


  print(map + map2 + map_se + map_se2 + plot_layout(design = c("
                                                               12
                                                               34"),
                                                    guides = "collect"))
  
  print(map + map2 + t_plot1 + t_plot2 + plot_layout(design = c("
                                                               12
                                                               34")))
  
  dev.off()
  