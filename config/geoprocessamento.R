repair_shp <- function(shp) {
  if (!(shp %>% is("Spatial"))){
    stop("Arquivo da área de estudo inválido!")
  }
  
  if (is.na(crs(shp))){
    crs(shp) <- CRS("+init=epsg:4326")
  } 
  if (class(shp) =="SpatialPolygonsDataFrame") {
    shp <- shp %>% 
      rgeos::gBuffer(byid=TRUE, width=0) 
  } 
  return(shp)
}

area_map <- function(shp, title="", crs_subtitle=T, lat="lat", long="long", group="group", colour="black", fill=NA) {
  number_format <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)
  
  map_data <- NULL
  ifelse(is(shp, "Spatial"), map_data <- fortify(shp), map_data <- shp)
  
  map_tmp <- ggplot(
      data =  map_data,
      aes_string(
        x = long, 
        y = lat, 
        group = group
      )
    ) + 
    scale_x_continuous(labels = number_format) +
    scale_y_continuous(labels = number_format) +
    coord_equal()

  if (title!=""){
    map_tmp <- map_tmp + 
      labs(title=title)
  } 
  
  if (crs_subtitle && is(shp, "Spatial")){
    map_tmp <- map_tmp + 
      labs(subtitle = paste0(crs(shp)))
  }
  
  if (is.na(fill)){
    map_tmp <- map_tmp + 
      geom_polygon(colour = colour, fill = NA)
  } else {
    map_tmp <- map_tmp + 
      geom_polygon(colour = NA, aes_string(fill = fill))
  }

  return(map_tmp)  
}

greater_than <- function(shp_area, inf_limit=0) {
  if (!(shp_area %>% is("Spatial"))){
    stop("Arquivo da área de estudo inválido!")
  }
  if (inf_limit<=0){
    stop("Limite de área de inválido!")
  }
  
  if (class(shp_area) =="SpatialPolygonsDataFrame") {
    shp_tmp <- shp_area %>% 
      disaggregate()
    
    kept_areas <- shp_tmp[raster::area(shp_tmp)>=inf_limit, ]

    if (length(kept_areas@polygons)>0) {
      return(kept_areas)
    }
  }
  return(shp_area)
}


#shp = shape_area_estudo 
#cell_width = largura_celula 
#cell_height = altura_celula
#var_names = nomes_variaveis_grid


make_grid <- function(shp, cell_width=0, cell_height=0, var_names=NULL, centroid=T){
  if (!(shp %>% is("Spatial"))){
    stop("Arquivo da área de estudo inválido!")
  }
  if (cell_width<=0 || cell_height<=0){
    stop("Tamanho de célula inválido!")
  }
  
  shp@data <- shp@data %>%
    rename_all(tolower)

  if (!is.null(var_names)){
      shp@data <- shp@data %>%
        dplyr::select(var_names %>% tolower() %>% all_of()) 
  } 
  
  shp_area_bkp <- shp
  if (is(shp, "SpatialPolygonsDataFrame")){
    #shp_area_bkp <- shp_area
    shp <- shp %>% 
      #raster::buffer(
      #  (c(cell_width, cell_height) %>% min()), 
      #  dissolve=T
      #) %>%
      rgeos::gBuffer(width=0, byid=T) %>%
      as("SpatialPolygonsDataFrame")
  } 
  
  shp@data <- shp@data %>%
    rowid_to_column("cell_id")
  
  shp_tmp_file <- tempfile() %>% paste0(".shp")
  shp %>%
    writeOGR(
      dsn = shp_tmp_file, 
      layer=".", 
      driver="ESRI Shapefile",
      overwrite=T
    )
  
  raster_tmp_file <- tempfile() %>% paste0(".tif")
  raster_area <- shp %>%
    raster::raster(vals=0) %>% 
    writeRaster(
      raster_tmp_file,
      format = "GTiff",
      bylayer = F, 
      options = c("dstnodata =-9999.0"),
      overwrite=T
    )

  
  #shp_grid <- gdal_rasterize(
  #  shp_tmp_file,
  #  raster_tmp_file,
  #  burn = 0,
  #  at = T, 
  #  co = c("BIGTIFF=YES"), 
  #  a_nodata = "-9999.0",    
  #  tr = c(cell_width, cell_height),
  #  tap= T,
  #  ot = 'Float32',
  #  output_Raster = T,
  #  te = shp %>% bbox() %>% as("vector"),
  #  verbose = T
  #) %>% 
  #  rasterToPolygons() 
  # https://gis.stackexchange.com/questions/166753/fastest-way-to-convert-big-raster-to-polyline-using-r-or-python/313550#313550
  # https://gis.stackexchange.com/questions/192771/how-to-speed-up-raster-to-polygon-conversion-in-r/357792#357792
  
  shp_grid <- terra::rasterize(
      shp[,1],
      raster_area,
      touches = T, 
      field=shp[,1]@data, na.rm=F)
  
  values(shp_grid) <- seq(1:ncell(shp_grid))
  
  shp_grid <- rasterToPolygons(shp_grid, na.rm=F) 
  
  crs(shp_grid) <- crs(shp)

  
  grid_cells <- shp_grid %>%
    gIntersects(
      shp_area_bkp, 
      byid = TRUE, 
      prepared = T, 
      returnDense = F
    ) %>%
    compact()
  
  ###
  shp_grid <- shp_grid[names(grid_cells) %>% as.integer(),] %>% 
    spChFIDs(as.character(1:length(grid_cells)))
  ###
  
  names(shp_grid) <- "cell_id"
  
  shp <- shp_grid
  
  if (centroid){
    centroids <- shp %>% 
      gCentroid(byid=TRUE)
    
    shp@data <- shp@data %>%
      bind_cols(as.data.frame(centroids@coords)) %>%
      rename(x_centroid = x, y_centroid = y)
  }
  return(shp)
}



add_raster <- function(shp, raster_folder=NULL, var_names=NULL, scenario=NULL){
  if (!(shp %>% is("Spatial"))){
    stop("Arquivo da área de estudo inválido!")
  }
  if (is.null(raster_folder)){
    stop("Pasta de rasters inválida!")
  }
  if (is.null(var_names)){
    var_names <- raster_folder %>%
      list.files()
  } 
  if(is.null(scenario)){ # if scenario == NULL, we are using current data
    l <- raster_folder %>% list.files(full.names = T)
    var_names2 <- paste0(var_names, '\\D')
    raster_list <- l[grep(paste(var_names2,collapse="|"),l)]
      
    if (length(raster_list)!=length(var_names)){
      stop("Não é possível determinar com exatidão quais rasters serão usados a partir da lista de nomes informada!")
    }
    
    raster_stack <- raster_list %>%
      stack() %>%
      tryCatch(error = function(e) e)
    
  } else { # if scenario != NULL, we are using future data
    l <- raster_folder %>% list.files(full.names = T)
    raster_stack <- l[grep(scenario,l)] %>%
      stack() %>%
      tryCatch(error = function(e) e)
  }
  
  # if tryCatch returns error, than we need to work on raster data:
  # 1. arrumar todos os rasters (extensão/resolução/...)
  if(any(class(raster_stack)=='error')){
    stop("Rasters diferentes não podem ser juntos em um stack.")
  }
  
  if(!is.null(scenario)){
    names(raster_stack) <- paste0('bio_', 1:19)
    raster_stack <- raster_stack[[unlist(var_names)]]
  } else {
    raster_stack <- raster_stack[[str_sort(names(raster_stack), numeric=T)]]
    names(raster_stack) <- var_names
  }
  
  # 2. Obter bbox do shape
  # 3. crop mask stack
  raster_stack <- raster_stack %>%
    crop(shp) %>%
    mask(shp) %>%
    stack()
  
  v <- list()
  for(i in shp$cell_id){
    sc <- subset(shp, cell_id==i)
    vals <- as.vector(terra::extract(raster_stack, sc, fun=mean, na.rm=T))
    v[[i]] <- vals
  }
  
  v <- v %>% 
    do.call(rbind,.) %>% 
    as.data.frame()
  
  # 7. rasteriza shape (dentro do stars)
  st_rasterize(st_as_sf(cbind(shp, v[1:28,])))
  
  shp_raster <- shp %>%
    cbind(v) %>%
    st_as_sf() %>%
    st_rasterize()
  
  names(shp_raster)[-(1:3)] <- names(raster_stack)
  
  return(shp_raster)
}


