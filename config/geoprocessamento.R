
areas_greater_than <- function(shp_area, inf_limit){
  if (class(shp_area) =="SpatialPolygonsDataFrame"){
    shp_tmp <- shp_area %>% 
      disaggregate()
    
    kept_areas <- shp_tmp[raster::area(shp_tmp)>=inf_limit, ]

    if (length(kept_areas@polygons)>0){
      shp_area <- kept_areas
    }
  }
  return(shp_area)
}

make_grid <- function(shp_area, cell_width, cell_height, var_names, withCentroid=T){
  shp_area@data <- shp_area@data %>%
    rename_all(tolower) %>%
    select(var_names %>% tolower() %>% all_of()) 
  
  # Funciona, mas é lento, não usado po enquanto, colocado aqui apenas para lembrar como faz
  .make_grid_lines_rasterize <- function(shp_area, cell_width, cell_height){
    .cell_count <- 0
    .cell_list_df <- list()
    .shp_data <- shp_area@data %>%
      rowid_to_column("id_celula")
    
    .lineshp_to_raster <- function(x, ...){
      .cell_count <<- .cell_count + 1
      if (!anyNA(x)){
        x <- .shp_data %>% 
          filter(id_celula %in% x) %>%
          select(-id_celula)
        .cell_list_df[[.cell_count]] <<- x %>%
          colMeans()
        return(.cell_count)
      } else{
        return(NA)
      }
    }
    
    .raster_mask <- shp_area %>% 
      extent() %>% 
      raster(resolution=c(cell_width, cell_height), crs = crs(shp_area)) %>%
      extend(c(1,1))
    
    beginCluster(n = ifelse(detectCores()-1>8, 8,detectCores()-1))
    shp_grid <- shp_area %>% 
      rasterize(.raster_mask, fun=.lineshp_to_raster) %>% 
      rasterToPolygons()
    endCluster()
    
    shp_grid@data <- .cell_list_df %>% 
      map_dfr(~ bind_rows(.)) 
    
    return(shp_grid)
  }
  
  # Extremamente rápido, mas não funciona para poligonos muito grandes, como todo o globo terrestre.
  # Funciona muito bem para linhas
  .make_grid_intersect <- function(shp_area, cell_width, cell_height){
    shp_grid <- shp_area %>% 
      extent() %>% 
      raster(
        resolution=c(cell_width, cell_height), 
        crs = crs(shp_area)
      ) %>%
      extend(c(1,1)) %>%
      rasterToPolygons()
    
    grid_cells <- shp_grid %>%
      gIntersects(
        shp_area, 
        byid = TRUE, 
        prepared=T, 
        returnDense = F
      ) %>%
      compact()
    
    shp_grid <- shp_grid[names(grid_cells) %>% as.integer(),] %>% 
      spChFIDs(as.character(1:length(grid_cells)))
    
    shp_grid@data <- grid_cells %>%
      map_dfr(as_tibble, .id="id_celula") %>% 
      rename(line_id=value) %>% 
      left_join(
        shp_area@data %>% 
          rowid_to_column("line_id")
      ) %>%
      select(-line_id) %>%
      mutate(id_celula=as.integer(id_celula)) %>%
      group_by(id_celula) %>% 
      summarise_all(mean) %>%
      ungroup() %>%
      select(-id_celula) %>%
      rowid_to_column("id_celula")
    
    return(shp_grid)
  }
  
  #funciona muito bem para poligonos inclusive grandes
  .make_grid_pol_rasterize <- function(shp_area, cell_width, cell_height){
    shp_area@data <- shp_area@data %>%
      rowid_to_column("id_celula")
    
    .raster_mask <- shp_area %>% 
      extent() %>% 
      raster(resolution=c(cell_width/10, cell_height/10), crs = crs(shp_area)) %>%
      extend(c(1,1))
    
    shp_grid <- shp_area %>% 
      fasterizeDT(.raster_mask, field = "id_celula") %>% 
      aggregate(10, fun=median) %>%
      rasterToPolygons() 
    
    shp_grid@data <- shp_grid@data %>% 
      rename(id_celula=layer) %>% 
      left_join(shp_area@data)
    
    shp_grid@data$id_celula <- 1:nrow(shp_grid@data)
    
    shp_grid <- shp_grid %>% 
      spChFIDs(as.character(1:nrow(shp_grid@data)))
    
    return(shp_grid)
  }
  
  if (is(shp_area, "SpatialLines")){
    shp_grid <- .make_grid_intersect(shp_area, cell_width, cell_height)
  } else {
    shp_grid <- .make_grid_pol_rasterize(shp_area, cell_width, cell_height)
  }
  
  if (withCentroid){
    centroids <- shp_grid %>% 
      gCentroid(byid=TRUE)
    
    shp_grid@data <- shp_grid@data %>%
      bind_cols(as.data.frame(centroids@coords))
  }
  return (shp_grid)
}

add_raster <- function(shp, raster_filenames, var_names, cell_width, cell_height){
  brick_reescaled <- brick() 
  for (i in 1:length(raster_filenames)){
    raster_file <- raster_filenames[[i]] %>%
      raster()
    
    raster_temp_file <- tempfile(fileext = '.tif')
    
    if (is.na(crs(raster_file))){
      crs_temp <- CRS("+init=epsg:4326")
    } else {
      crs_temp <- crs(raster_file)
    }
    
    raster_rescaled <- gdalwarp(
      raster_filenames[[i]],
      raster_temp_file,
      s_srs = crs_temp,
      t_srs = crs(shp),
      te = shp %>% bbox() %>% as.vector(),
      r = 'lanczos',
      tr = c(cell_width, cell_height),
      ot = 'Float32',
      multi = T,
      output_Raster = T,
      overwrite = T
    )
    
    raster_rescaled[is.na(raster_rescaled[])] <- 0
    brick_reescaled  <- brick_reescaled %>% addLayer(raster_rescaled)
  }
  
  brick_reescaled <- brick_reescaled %>% 
    mask(shp)
  
  new_shp <- brick_reescaled %>% 
    rasterToPolygons()
  
  names(new_shp) <- var_names
  file_delete(raster_temp_file)
  
  shp@data <- shp@data %>%
    bind_cols(new_shp@data) %>%
    select(c("id_celula", var_names) %>% all_of()) %>%
    mutate(id_celula=1:nrow(shp@data)) %>% 
    as.data.frame()
    
  shp <- shp %>% 
    spChFIDs(as.character(1:nrow(shp@data)))
  
  return(shp)
}


repair_shapefile <- function(shp){
  if (class(shp) =="SpatialPolygonsDataFrame"){
    shp %>% 
      gBuffer(byid=TRUE, width=0) %>%
      return()
  } else{
    return(shp)
  }
}

reproject_shapefile <- function(shp, proj){
  if (!(shp %>% crs() %>% compareCRS(CRS(proj)))){
    shp %>% 
      spTransform(CRS(proj)) %>%
      return()
  } else{
    return(shp)
  }
}


map_of_area <- function(shp, title="", crs_subtitle=T, lat="lat", long="long", group="group", colour= "black", fill=NA){
  if (is(shp, "Spatial")){
    map_tmp <- ggplot(
        data =  fortify(shp),
        aes_string(
          x = long, 
          y = lat, 
          group = group
        )
    )
  } else {
    map_tmp <- ggplot(
      data =  shp,
      aes_string(
        x = long, 
        y = lat, 
        group = group
      )
    )
  }
  
  map_tmp <- map_tmp +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma) +
    coord_equal()
  
  if (title!=""){
    map_tmp <- map_tmp + labs(title=title)
  } 
  
  if (crs_subtitle && is(shp, "Spatial")){
    map_tmp <- map_tmp + labs(subtitle = paste0(crs(shp)))
  }
  

  if (is.na(fill)){
    map_tmp <- map_tmp + geom_polygon(colour = colour, fill = NA)
  } else {
    map_tmp <- map_tmp + geom_polygon(colour = NA, aes_string(fill = fill))
  }

  return(map_tmp)  
}

