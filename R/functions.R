# Calculate KDE at multiple levels -------------------------
gen_kde <- function(points, strata){
  # if(!is_null(rm_rows)) points <- points[-rm_rows,]
  sr <- points |> group_by({{strata}})
  s <- sr |> 
    group_split()
  birdid <- group_keys(sr)[[1]]
  # site <- gr_tri_loc$Location[gr_tri_loc$BirdID==birdid]
  # browser()
  names(s) <- birdid
  # bw <- map(s, ~)
  # browser()
  r_kde <- map(1:length(s),#bw, 
               ~{
                 tmp <- (
                   eks::st_kde(s[[.x]],
                               H=ks::Hlscv(distinct(as_tibble(st_coordinates(s[[.x]])), X, Y),
                                           optim.fun = 'nlm') ))
                 tmp$sf$BirdID <- names(s)[[.x]]
                 tmp$sf$Area <- st_area(tmp$sf$geometry)
                 tmp$sf$contlabel <- factor(tmp$sf$contlabel)
                 tmp
               }
               
  )
  # r <- map(r_kde, raster::raster)
  names(r_kde)<- names(s)
  return(r_kde)
  
}

jack_kde <- function(points_df){
  res <- map(1:nrow(points_df), ~{
    points <- points_df[-.x,]
    gen_kde(points, strata = BirdID)
  } )
  return(res)
}


# Pull out the 95% KDEs
pull95 <- function(x,list_, list_for_names){
  list_[[x]] |> 
    transpose() |> 
    pluck(names(list_for_names)[[x]]) |>  
    transpose() |> 
    pluck("sf") |> 
    bind_rows() |> 
    filter(contlabel==95) }


fgg <-  function(.x, df){
  df |> 
    filter(BirdID==.x) |> 
    left_join(band_ids, by = "BirdID") |> 
    ggplot()+ 
    geom_sf(alpha = 0.5,
            fill = NA)  +
    # labs(title = .x) +
    facet_wrap(Location~TagID)+
    theme_minimal()
}



get_overlap <- function(df, Location_, type_=NULL, contlabel_=NULL, byBird=F){
  if(isTRUE(byBird)){
    comp <- tibble(.x = c("Visual", "Visual", "Triangulation"),
                   .y = c("Triangulation", "Home Range", "Home Range") )
    comp$compare_area <- pmap(comp,
                              ~sum(st_area(st_intersection(df |> filter(type ==.x),
                                                           df |> filter(type==.y)))) ) |> list_c()
    comp$BirdID <- Location_
    comp$comp_lab <- glue::glue("{comp$.x}_{comp$.y}")
    return(comp)
  }else{
    d <- filter(df, Location==Location_, type==type_, contlabel==contlabel_)
    
    d$overlap <- 
      map(1:nrow(d), ~sum(st_area(st_intersection(d[.x,], d[-.x,])))) |> list_c()
    
    d}
  
}



gen_overlap_plot <- function(df, Location_, type_, contlabel_){
  d <- filter(df, Location==Location_, type==type_, contlabel==contlabel_)
  
  
  
  ggplot() +
    geom_sf(data = d,
            fill = 'white', colour = "black") + #, aes(colour = TagID)) +
    geom_sf(data =st_intersection(d) |> 
              filter(n.overlaps>1),fill = 'darkgrey')+#aes(fill = n.overlaps)) +
    ggthemes::theme_map(base_size = 14, base_family = "Roboto Condensed") +
    scale_colour_viridis_d()+
    labs(colour = "Bird", title = ifelse(Location_=="Panmure", glue::glue(" {type_}"), "")) + #{Location_}:
    ggsn::scalebar(st.size = 2,st.dist = 0.04,border.size = .5,
                   data = d, dist_unit = 'm', transform = F, dist = 100, location = 'bottomleft', )
  
}


individual_plots <- function(.x){
  d_id <- 
    all_poly |> 
    filter(contlabel  %in% c(5, 50, 95) & TagID==.x &
             ((  type == "Home Range" & contlabel == 95) |
                (type != "Home Range" & contlabel ==95))
    ) 
  plt <- 
    d_id |> 
    ggplot() + geom_sf(aes( colour = type),
                       linewidth =1,fill = NA) +
    theme_minimal() +
    labs(title = .x) +
    geom_sf(data = VisualPoints |> 
              st_as_sf() |> 
              left_join(band_ids, by = "BirdID") |> 
              filter(TagID == .x)) +
    geom_sf(data = TriangulationPoints |> 
              st_as_sf() |>  
              filter(TagID == .x ), shape =2) +
    ggthemes::theme_map(base_size = 14, base_family = "Roboto Condensed") +
    rcartocolor::scale_colour_carto_d("", palette = 4)+
    # labs(colour = "Bird", title = ifelse(Location_=="Panmure", glue::glue(" {type_}"), "")) + #{Location_}:
    ggsn::scalebar(st.size = 2,st.dist = 0.04,border.size = .5,
                   data = d_id, dist_unit = 'm', transform = F,
                   dist = 50, location = 'bottomleft' ) +
    theme(legend.position = 'right')
  ggsave(glue::glue("output/{str_replace(.x, '/', '_')}_indiv.jpg"), plt)
  plt
}






sgen_kde <- safely(gen_kde)

run_bootstrap <- function(points, strata, n, seed){
  withr::with_seed(seed, {
    xy <- slice_sample(points,n=n, by = {{strata}}, replace = T)
  })
  
  sgen_kde(xy, {{strata}})$result %>% {
    if(is.null(.)) {tibble(BirdID = NA) }
    else{
      transpose(.) |> 
        pluck('sf') |> 
        bind_rows() |> 
        filter(contlabel == 95) |> 
        st_drop_geometry() |> 
        dplyr::select(BirdID, Area) 
    }}
  
  
}