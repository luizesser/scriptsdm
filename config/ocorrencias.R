occurrences_to_shapefile <- function(o_file, sp_names, shp, proj){
  sp_names <- sp_names %>% 
    to_snake_case() %>% 
    abbreviate(minlength = 10)
  
  suppressMessages(
    sp_occurrences <- o_file %>%
      vroom() %>%
      clean_names() %>%
      remove_empty(c("rows", "cols")) %>%
      mutate(especie=especie %>% to_snake_case() %>% abbreviate(minlength = 10)) 
  )
  if (!is.na(sp_names) && !is.null(sp_names) && length(sp_names)>0)
    sp_occurrences <-  sp_occurrences %>%
      filter(especie  %in% sp_names)
  
  if (is.character(shp))
    shp <- shp %>% 
      readOGR(verbose = F)
  
  sp_occurrences$longitude <- sp_occurrences$longitude %>% 
    as.character() %>% 
    as.numeric()
  
  sp_occurrences$latitude <- sp_occurrences$latitude %>% 
    as.character() %>% 
    as.numeric()
  
  coordinates(sp_occurrences) <- ~longitude+latitude
  if (is.character(proj)){
    crs(sp_occurrences) <- CRS(proj)
  } else {
    crs(sp_occurrences) <- proj  
  }
  
  sp_occurrences %>% 
    spTransform(crs(shp)) %>%
    return()
}

occurrences_to_pa_shapefile <- function(shp_occ, shp_area, sp_names){
  sp_names <- sp_names %>% 
    to_snake_case() %>% 
    abbreviate(minlength = 10)
  
  pa_matrix <- shp_area@data[,F]
  
  for (spp in sp_names){
    pa_matrix <- pa_matrix %>% 
      bind_cols(
        shp_area %>%
          over(shp_occ[shp_occ@data$especie==spp, ])
      )
  }
  names(pa_matrix) <- sp_names
  
  pa_matrix[!is.na(pa_matrix)] <- 1
  pa_matrix[is.na(pa_matrix)] <- 0
  
  pa_matrix <- pa_matrix %>% 
    mutate_all(~ as.integer(.))
  
  rownames(pa_matrix) <- rownames(shp_area@data)
  
  grid_pa_matrix <- shp_area 
  grid_pa_matrix@data <-  pa_matrix
  
  return(grid_pa_matrix)
}

map_of_occurrences <- function(shp_occ, shp_area, title="", crs_subtitle=T){
  map_tmp <- ggplot(
    data =  fortify(shp_area),
    aes(
      x = long, 
      y = lat,
      group=group
    )
  )

  map_tmp <- map_tmp +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma) +
    coord_equal()
  
  if (title!=""){
    map_tmp <- map_tmp + 
      labs(title=title)
  } 
  
  if (crs_subtitle){
    map_tmp <- map_tmp + 
      labs(subtitle = paste0(crs(shp_area)))
  }

  map_tmp <- map_tmp + 
    geom_polygon(colour = "black", fill=NA)
  
  map_tmp <- map_tmp + 
    geom_point(data = fortify(shp_occ), aes(x = long, y = lat, group=especie, color=especie)) 

  return(map_tmp)  
}

map_of_pa <- function(shp_pa, shp_area, sp_names){
  sp_names <- sp_names %>% 
    to_snake_case() %>% 
    abbreviate(minlength = 10)
  
  df_temp <- shp_pa@data %>% 
    rownames_to_column("id") %>%
    pivot_longer(-id, names_to = "especies", values_to="presenca") 
  
  legend_title  <- ifelse(length(sp_names)==1, sp_names[[1]], "Presença/Ausência")
  
  ggplot(data = fortify(shp_area) %>% left_join(df_temp)) +
       aes(x = long, y = lat, group = group) +
       geom_polygon(aes(fill = as.factor(presenca))) +
       scale_fill_manual(legend_title, values=c("darkseagreen2", "tomato")) +
       scale_x_continuous(labels = comma) +
       scale_y_continuous(labels = comma) +
       coord_equal()
}