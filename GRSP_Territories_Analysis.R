library(dplyr)
library(ggplot2)
library(units)
library(remotes)
library(radiotrack)
library(raster)
library(sp)
library(ks)
library(sf)
library(tidyverse)
library(eks)
library(patchwork)
library(ggthemes)
source("R/functions.R")

# Setup ----------------------------
## Get triangulated points form telemetry data ####
GRSP_tri<-read.csv("data/TelemetryPoints.csv") |> 
  as_tibble() |> 
  mutate(dt= ymd_hms(DateTime_, tz = 'America/Toronto'),
         date = ymd(glue::glue("2016-{Month_}-{Day_}")),
         lat = Lat, lon = Long_) %>%
  # rowwise() |> 
  mutate(sr = suncalc::getSunlightTimes(data = ., keep = "sunrise", tz = 'America/Toronto')[["sunrise"]],
         t2sr = difftime( dt,sr, units = 'mins'))

## Join location to band and bird IDs ----------------------
gr_tri_loc <- dplyr::select(GRSP_tri, BirdID = BIRD, Location) |> 
  distinct()

band_ids <-read_csv("data/Bird_band_ids.csv") |> 
  mutate(BirdID = stringr::str_to_upper(str_remove(nickname, "Papa")),
         TagID = paste0(Right_upper, "/", Right_lower)) |> 
  full_join(gr_tri_loc)



## Calculate locations based on triangulation ----------------------
results<- triangulate(GRSP_tri, x = "Longitude", y = "Latitude",
                      bearings = "Bearing", group = "ID")






## Final layers for locations ---------------------------------
TriangulationPoints <- results |> 
  filter(Error == "") |> 
  st_as_sf(coords = c("X_Coordinate", "Y_Coordinate"), crs = 4326) |> 
  st_transform(3161) |> 
  mutate(BirdID=str_extract(ID, "^[A-Z]+(?=(_))")) |> 
  # st_read(
  # "GIS/TriangulationPoints_latlong.shp") |>
  left_join(band_ids)
VisualPoints <- st_read(
  "data/VisualPoints_latlong.shp") |> 
  mutate(dt= ymd_hms(DateTime_, tz = 'America/Toronto'),
         date = ymd(glue::glue("2016-{Month_}-{Day_}")),
         lat = Lat, lon = Long_) %>%
  # rowwise() |> 
  mutate(sr = suncalc::getSunlightTimes(data = ., keep = "sunrise", tz = 'America/Toronto')[["sunrise"]],
         t2sr = difftime( dt,sr, units = 'mins')) |> 
  st_transform(3161) |> left_join(band_ids)

sites_s <- 
  GRSP_tri |> 
  st_as_sf(coords=c("Longitude", "Latitude"),crs= 4326) |> 
  st_transform(3161) |> 
  bind_rows( VisualPoints) |> 
  st_as_sf() |> 
  dplyr::select(Location) |> distinct() |>
  st_transform(3161) |> 
  st_buffer(100) |> 
  st_transform(4326) |> 
  group_split(Location, .keep = T) 
sites <- sites_s %>%
  map(., st_bbox)
names(sites) <- c(sites_s[[1]]$Location[[1]], sites_s[[2]]$Location[[1]]) # Manually add site names





# Analyses -------------------------
## Run kde estimation for Triangulation, Visual, and all points by Bird ----
res <- gen_kde(TriangulationPoints , strata = BirdID)
resV <- gen_kde(VisualPoints, strata = BirdID)
res_home_range <- gen_kde(bind_rows(TriangulationPoints,VisualPoints ) |> 
                            dplyr::select(BirdID, Location)
                          , strata = BirdID)



## Jackknife estimation ---------------------
TriangulationPoints_sep <- split(x = TriangulationPoints, f = TriangulationPoints$BirdID, drop = FALSE)
VisualPoints_sep <- split(x = VisualPoints, f = VisualPoints$BirdID, drop = FALSE)

HR_dat <- bind_rows(TriangulationPoints,VisualPoints ) |> 
  dplyr::select(BirdID, Location)
HomeRange_sep <- split(x = HR_dat, f = HR_dat$BirdID, drop = FALSE)


jackknife_res_tri <- map(TriangulationPoints_sep, jack_kde)
jackknife_res_Vis <- map(VisualPoints_sep, jack_kde)
jackknife_res_HR <- map(HomeRange_sep, jack_kde)



res95jack_tri <- map(1:length(TriangulationPoints_sep),
                     pull95, list_ = jackknife_res_tri,
                     list_for_names = TriangulationPoints_sep) |> 
  list_rbind() |> st_as_sf()  |> mutate(type = "Triangulation")

res95jack_vis <- map(1:length(VisualPoints_sep),
                     pull95, list_= jackknife_res_Vis,
                     list_for_names = VisualPoints_sep) |> 
  list_rbind() |> st_as_sf() |> mutate(type = "Territory")

res95jack_hr <- map(1:length(HomeRange_sep),
                     pull95, list_= jackknife_res_HR,
                     list_for_names = HomeRange_sep) |> 
  list_rbind() |> st_as_sf() |> mutate(type = "Home Range")

mapJack_tri <- map(unique(res95jack_tri$BirdID),fgg, df = res95jack_tri)
mapJack_vis <- map(unique(res95jack_vis$BirdID),fgg, df = res95jack_vis)
mapJack_hr <- map(unique(res95jack_hr$BirdID),fgg, df = res95jack_hr)

pdf("output/Jackknife_Visual_plots.pdf")
mapJack_vis
dev.off()

pdf("output/Jackknife_Tri_plots.pdf")
mapJack_tri
dev.off()


pdf("output/Jackknife_HR_plots.pdf")
mapJack_hr
dev.off()






## Get 95% KDE polygons #### ---------------------



poly_vis <- transpose(resV) |> pluck("sf") |> bind_rows() |> 
  mutate(type = "Territory")
poly_tri <- transpose(res) |> pluck("sf") |> bind_rows() |> 
  mutate(type = "Triangulation")

poly_home_range <- transpose(res_home_range) |> pluck("sf") |> bind_rows() |> 
  mutate(type = "Home Range")


all_poly <- bind_rows(
  list(poly_vis, poly_tri, poly_home_range)
) |> left_join(band_ids,by = join_by(BirdID))

## Calculate area by bird -----------
comparison_by_bird <- 
  st_drop_geometry(all_poly) |> 
  dplyr::select(contlabel, BirdID,TagID, Area, type) |> 
  pivot_wider(names_from = type, values_from = Area) |> 
  filter(contlabel==95) |> 
  mutate(across(4:6, ~set_units(.x,"ha")))




### Plot outlining how Area changes with KDE level -------------
ggplot(all_poly, aes(as.numeric(as.character(contlabel)),
                     set_units(Area, "ha") , colour = type)) +
  geom_line() + facet_wrap(Location~TagID) +
  theme_minimal(base_size = 14, base_family = "Roboto Condensed" ) +
  labs(x = "KDE level", y = "Area", colour = "") +
  rcartocolor::scale_colour_carto_d()


ggsave("output/Change_Area_with_KDE.png", width = 12, height = 12, units = 'in', dpi = 600)

ggplot(all_poly |> 
         filter(contlabel == 95),
       aes(TagID,set_units(Area, "ha") , colour = type)) +
  
  stat_summary(data = bind_rows(res95jack_tri,
                                res95jack_vis, res95jack_hr) |> 
                 left_join(band_ids),
               # fun.data = 'mean_cl_boot', 
               fun.min = min,
               fun.max = max,
               geom = "linerange",
               position = position_dodge(width = 0.33))+
  geom_point(position = position_dodge(width = 0.33)) +
  facet_wrap(~Location, scales = 'free_x') +
  labs(x = "Bird tags", y = "Area with jackknife error",
       colour = "Estimation method") +
  rcartocolor::scale_colour_carto_d() +
  ggthemes::theme_clean(base_size = 14, base_family = "Roboto Condensed" ) + 
  guides(x = guide_axis(check.overlap = T, n.dodge = 2))

ggsave("output/Area_with_error.png", width = 12, height = 6, units = 'in', dpi = 600)




overlaps <- 
  all_poly |>
  dplyr::distinct(Location, type, contlabel) %>% 
  rename_with(~paste0(., "_"), everything()) |> 
  pmap(get_overlap, df = all_poly) |> 
  list_rbind() 


overlaps_by_bird <- 
  all_poly |>
  filter(contlabel==95) |> 
  mutate(Location_ = BirdID) |> 
  dplyr::select(type, BirdID) |> 
  nest_by(BirdID) |> 
  rowwise() |> 
  mutate(o= list(get_overlap( df = data,Location_=.data$BirdID,byBird=T))) |> 
  dplyr::select(o) |> unnest(o) |> 
  mutate(area_ha = set_units(compare_area, "ha")) |> 
  left_join(band_ids)

write_csv(
  dplyr::select(overlaps_by_bird,BirdID,TagID, Location, comp_lab, area_ha ),
  "output/Within_bird_overlap.csv"
)



  

## Table 1 ------------


Table1 <- 
  comparison_by_bird |> 
  filter(contlabel==95) |> 
  left_join(VisualPoints |> st_drop_geometry() |> summarise(nVisual=n(), .by = BirdID)) |> 
  left_join(TriangulationPoints |> st_drop_geometry() |> summarise(ntri=n(), .by = BirdID)) |> 
  mutate(across(4:6, ~set_units(.x,"ha"))) |> 
  mutate("Home Range points" = nVisual+ntri ) |> 
  left_join(band_ids) |> 
  left_join(overlaps_by_bird |> filter(comp_lab == "Visual_Home Range") |> 
              dplyr::select(BirdID, overlap_ha =area_ha)) |> 
  mutate("KDE overlap %" = overlap_ha/`Home Range`) |> 
  dplyr::select("Study Site" = Location,BirdID,
                TagID, "Territory KDE (ha)" =Territory,
                "Territory points" = nVisual,
                "Telemetry KDE (ha)"= Triangulation,
                "Telemetry points" = ntri,
                "Home Range KDE (ha)" = `Home Range`,
                "Home Range points", 
                "KDE overlap %" 
  ) 

Table1_dropped <- 
  Table1

Table1_dropped$`Territory KDE (ha)`[Table1_dropped$`Territory points`<10] <- NA
Table1_dropped$`Home Range KDE (ha)`[Table1_dropped$`Home Range points`<15] <- NA


gt::gt(Table1)

summarize(Table1, across(where(is.numeric), mean) )
summarize(Table1_dropped, across(where(is.numeric), ~mean(.x, na.rm=T) ))


summarize(Table1_dropped, across(where(is.numeric), ~sd(.x, na.rm=T)/sqrt(n()) )) |> 
  glimpse()







n_obs <- 
VisualPoints |> st_drop_geometry() |> summarise(n=n(), .by = BirdID) |> 
  mutate(type = "Visual") |> 
  bind_rows(
    bind_rows(VisualPoints,
    TriangulationPoints) |> st_drop_geometry() |> summarise(n=n(), .by = BirdID) |> 
      mutate(type = "Home Range") )








write_csv(Table1
  ,
  "output/Table1.csv")


## Overlaps ----------------


overlapplots <- 
  all_poly |>
  filter(contlabel == 95) |> 
  mutate(type = str_replace(type, "Territory", "Visual"))  |> 
  dplyr::distinct(Location, type, contlabel) %>% 
  rename_with(~paste0(., "_"), everything()) |> 
  pmap(gen_overlap_plot, df = all_poly) 

figure3 <- 
(wrap_elements(grid::textGrob("Panmure", check.overlap = F)) + 
   overlapplots[[1]] + 
    overlapplots[[3]] +overlapplots[[5]]  + 
    plot_layout(nrow =1, widths = c( 0.12, 0.3,0.3,0.3))) /
  (wrap_elements(grid::textGrob("Burnt\nLands", check.overlap = F) )+
     overlapplots[[2]] + overlapplots[[4]] +overlapplots[[6]] + 
     plot_layout(nrow =1,  widths = c( 0.12, 0.3,0.3,0.3)) )

ggsave("output/figure3.png",plot = figure3, width = 2000,
       height = 1000, units = 'px',
       dpi = 300)




overlaps |> 
  dplyr::select(Location, TagID, contlabel, Area,type, Overlap_area = overlap) |> 
  mutate(across(c(Area, Overlap_area), ~set_units(.x,"ha"))) |> 
  filter(contlabel == 95) |> write_csv("output/overlap_by_bird_ha.csv")



indiv_plts <- map(unique(band_ids$TagID[-13]), individual_plots)

individual_plots("Black/Orange")
file.copy("output/Black_Orange_indiv.png", "output/figure2.png", overwrite = T)

site_plots <- map(c("Burnt Lands", "Panmure"), ~{
  fd_ <- 
    all_poly |> 
    filter(contlabel  %in% c(95) & 
             Location == str_remove(.x, "\\s") &
             type !="Home Range") 
  fd_ |> 
    ggplot() +
    geom_sf(fill = NA,
            aes(linetype = type, colour = TagID),
            linewidth = .8) +
    ggthemes::theme_map(base_size = 12, base_family = "Roboto Condensed") +
    labs(linetype = "Detection\nmethod",
         colour = "Bands", title = .x) +
    theme(legend.position = 'right') +
    rcartocolor::scale_colour_carto_d()+ 
    ggsn::scalebar(st.size = 2,st.dist = 0.04,
                   border.size = .5,
                   data = fd_, dist_unit = 'm', 
                   transform = F,
                   dist = 100, location = 'bottomleft' ) 
})
ggsave(filename ="output/BurntLands.jpg",plot = site_plots[[1]] , width = 4, height = 6, dpi = 1200)
ggsave(filename ="output/Panmure.jpg",plot = site_plots[[2]] , width = 4, height = 6, dpi = 1200)


comp <- 
  comparison_by_bird |>
  filter(contlabel==95) |>
  left_join(VisualPoints |> st_drop_geometry() |> summarise(nVisual=n(), .by = BirdID)) |>
  left_join(TriangulationPoints |> st_drop_geometry() |> summarise(ntri=n(), .by = BirdID)) |>
  mutate(across(4:6, ~set_units(.x,"ha"))
  ) |> 
  pivot_longer(names_to = "Type", values_to = "ha", cols = c(Territory, Triangulation, `Home Range`),
               values_transform = units::drop_units) |> 
  left_join(band_ids)




# Run bootstraps ------------



res_home_range <- gen_kde(bind_rows(TriangulationPoints,VisualPoints ) |> 
                            dplyr::select(BirdID, Location), 
                          strata = BirdID)

bootstraps <- 
  expand_grid(n_=2:50, iter=1:1000) |> 
  rowwise() |> 
  mutate(res = list(run_bootstrap(
    bind_rows(TriangulationPoints,VisualPoints ) |>
      dplyr::select(BirdID, Location),
    strata = BirdID, n = n_, seed = iter) ) ) |>
  ungroup() |> unnest(res)

write_rds(bootstraps, "output/Boostraps.rds")


d <- read_rds("output/Boostraps.rds")
glimpse(d)
# d$Area |> n_distinct()

d |> filter(!is.na(Area) & n_%%5==0 ) |> 
  mutate(across(Area,  ~units::set_units(.x, 'ha'))) |> 
  ggplot(aes(n_,Area, group = BirdID)) +
  # stat_summary(
  #              fun.data = 'mean_sdl', geom='ribbon', alpha = 0.4) + 
  stat_summary(
    fun.min = function(z) { quantile(z,0.025) },
    fun.max = function(z) { quantile(z,0.975) },
    fun = median,
    geom='ribbon', alpha = 0.4) +
  stat_summary(fun = 'mean', geom='line') +
  labs(x = "N observations per bird",
       y = expression(paste("Mean home range (mean"%+-%"95%CI)") ))+
  theme_light() + 
  geom_hline(data=Table1 , aes(yintercept = `Home Range KDE (ha)`), linetype =2) +
  geom_vline(data=Table1 , aes(xintercept = `Home Range points`), linetype =2) +
  facet_wrap(~BirdID, scales = 'free')



r <- 
  tibble(n_ =
           bind_rows(TriangulationPoints,VisualPoints ) |> 
           st_drop_geometry() |> 
           count(BirdID) |> 
           pull(n) |> mean(),
         Area = 
           poly_home_range |> 
           filter(contlabel==95) |> 
           st_drop_geometry() |> 
           pull(Area) |> 
           mean() |> 
           units::set_units('ha'))

bootstraps_plot <- 
d |> filter(!is.na(Area)) |> 
  mutate(across(Area,  ~units::set_units(.x, 'ha'))) |> 
  ggplot(aes(n_,Area)) +
  stat_summary(
    fun.data = 'mean_se',
    geom='ribbon', alpha = 0.4) +
  stat_summary(fun = 'mean', geom='line') +
  labs(x = "N observations per bird",
       y = expression(paste("Mean home range (mean"%+-%"SE)") ))+
  ggthemes::theme_clean(base_family = 'arial') +
  # geom_hline(yintercept = r$Area) +
  geom_vline(xintercept = r$n_)
ggsave("output/bootstraps_home_range.png", dpi = 300,units = 'px',
       width = 800, height = 800, plot = bootstraps_plot)


job::job({
  runs_terr <- expand_grid(n_=c(seq(2,19, by = 2), seq(20, 50, by = 5)), iter=1:1000) 
  
  library(furrr)
  plan(multisession(workers = 16))
  boot_res <- future_map(1:nrow(runs_terr),
                       ~run_bootstrap(
                         VisualPoints  |>
                           dplyr::select(BirdID, Location),
                         strata = BirdID, n = runs_terr$n_[[.x]],
                         seed = runs_terr$iter[[.x]]))

plan(sequential)
file.create("C:/Documents and Settings/hoped/OneDrive - EC-EC/MS_collaborations/GRSP/Bootstraps2_complete.txt")
})
rt <- bind_rows(boot_res)
runs_terr$res <- boot_res
write_rds(unnest(runs_terr,res ) |> 
            filter(!is.na(Area)), "output/Bootstraps_terr.rds")

d_terr<- read_rds("output/Bootstraps_terr.rds")

d_terr |> filter(!is.na(Area) & n_%%5==0 ) |> 
  mutate(across(Area,  ~units::set_units(.x, 'ha'))) |> 
  ggplot(aes(n_,Area, group = BirdID)) +
  # stat_summary(
  #              fun.data = 'mean_sdl', geom='ribbon', alpha = 0.4) + 
  stat_summary(
    fun.min = function(z) { quantile(z,0.025) },
    fun.max = function(z) { quantile(z,0.975) },
    fun = median,
    geom='ribbon', alpha = 0.4) +
  stat_summary(fun = 'mean', geom='line') +
  labs(x = "N observations per bird",
       y = expression(paste("Mean territory size (mean"%+-%"95%CI)") ))+
  theme_light() + 
  # geom_hline(data=Table1 , aes(yintercept = `Home Range KDE (ha)`), linetype =2) +
  # geom_vline(data=Table1 , aes(xintercept = `Home Range points`), linetype =2) +
  facet_wrap(~BirdID, scales = 'free')


r_terr <- 
  tibble(n_ =
           bind_rows(VisualPoints ) |> 
           st_drop_geometry() |> 
           count(BirdID) |> 
           filter(n>=10) |> 
           pull(n) |> mean(),
         Area = 
           poly_vis |> 
           filter(contlabel==95) |> 
           st_drop_geometry() |> 
           pull(Area) |> 
           mean() |> 
           units::set_units('ha'))
  

d_terr |> filter(!is.na(Area)) |> 
  mutate(across(Area,  ~units::set_units(.x, 'ha'))) |> 
  ggplot(aes(n_,Area)) +
  # stat_summary(
  #   fun.min = function(z) { quantile(z,0.025) },
  #   fun.max = function(z) { quantile(z,0.975) },
  #   fun = median,
  #   geom='ribbon', alpha = 0.4) +
  stat_summary(
    fun.data = 'mean_se',
    geom='ribbon', alpha = 0.4) +
  stat_summary(fun = 'mean', geom='line') +
  labs(x = "N observations per bird",
       y = expression(paste("Mean territory size (mean"%+-%"SE)") ))+
  theme_light() +
  # geom_hline(yintercept = r$Area) +
  geom_vline(xintercept = 10, linetype =2)+
  geom_vline(xintercept = r_terr$n_)


