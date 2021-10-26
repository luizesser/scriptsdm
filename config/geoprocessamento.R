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
  
  shp_grid <- gdal_rasterize(
    shp_tmp_file,
    raster_tmp_file,
    burn = 0,
    at = T, 
    co = c("BIGTIFF=YES"), 
    a_nodata = "-9999.0",    
    tr = c(cell_width, cell_height),
    tap= T,
    ot = 'Float32',
    output_Raster = T,
    te = shp %>% bbox() %>% as("vector"),
    verbose = T
  ) %>% 
    rasterToPolygons() 
  # https://gis.stackexchange.com/questions/166753/fastest-way-to-convert-big-raster-to-polyline-using-r-or-python/313550#313550
  # https://gis.stackexchange.com/questions/192771/how-to-speed-up-raster-to-polygon-conversion-in-r/357792#357792
  
  crs(shp_grid) <- crs(shp)
  
  shp_grid@data <- shp_grid@data %>%
    rowid_to_column("cell_id")
  
  grid_cells <- shp_grid %>%
    gIntersects(
      shp_area_bkp, 
      byid = TRUE, 
      prepared = T, 
      returnDense = F
    ) %>%
    compact()
  
  shp_grid <- shp_grid[names(grid_cells) %>% as.integer(),] %>% 
    spChFIDs(as.character(1:length(grid_cells)))
  
  shp_grid@data <- grid_cells %>%
    map_dfr(as_tibble, .id="cell_id") %>% 
    rename(line_id=value) %>% 
    left_join(
      shp_area_bkp@data %>% 
        rowid_to_column("line_id")
    ) %>%
    select(-line_id) %>%
    mutate(cell_id=as.integer(cell_id)) %>%
    group_by(cell_id) %>% 
    summarise_all(~ ifelse(is.numeric(.), mean(.), .), na.rm = TRUE) %>%
    ungroup() %>%
    select(-cell_id) %>%
    rowid_to_column("cell_id")

  
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

add_raster <- function(shp, raster_folder=NULL, cell_width=0, cell_height=0, var_names=NULL){
  if (!(shp %>% is("Spatial"))){
    stop("Arquivo da área de estudo inválido!")
  }
  if (is.null(raster_folder)){
    stop("Pasta de rasters inválida!")
  }
  if (cell_width<=0 || cell_height<=0){
    stop("Tamanho de célula inválido!")
  }
  
  if (is.null(var_names)){
    var_names <- raster_folder %>%
      list.files()
  } 
  
  raster_list <- raster_folder %>%
    list.files(full.names = T) %>% 
    keep(~ .x %>% str_detect(fixed(var_names, ignore_case = T)) %>% any())
    
  if (length(raster_list)!=length(var_names)){
    stop("Não é possível determinar com exatidão quais rasters serão usados a partir da lista de nomes informada!")
  }
  
  raster_stack <- raster_list %>%
    stack()

  raster_tmp_file <- tempfile() %>% paste0(".tif")
  raster_stack %>% 
    writeRaster(
      raster_tmp_file,
      format = "GTiff",
      bylayer = F, #bylayer = T,
      #suffix = lista_rasters_bio %>% names(),
      options = c("dstnodata =-9999.0"),
      overwrite=T,
    )
  
  shp_countour_file <- tempfile() %>% paste0(".shp")
  shp %>%
    raster::aggregate(dissolve=T) %>%
    #rgeos::gBuffer(width=-(min(c(cell_width, cell_height)))/10, capStyle = "SQUARE", joinStyle = "BEVEL") %>%
    #rgeos::gUnionCascaded() %>%
    as("SpatialPolygonsDataFrame") %>% 
    writeOGR(
      dsn = shp_countour_file, 
      layer=".", 
      driver="ESRI Shapefile",
      overwrite=T
    )
  
  shp_area_file <- tempfile() %>% paste0(".shp")
  shp %>%
    writeOGR(
      dsn= shp_area_file, 
      layer=".", 
      driver="ESRI Shapefile",
      overwrite=T
    )
  
  shp_grid_file <- tempfile() %>% paste0(".shp")
  shp %>%
    as("SpatialPolygonsDataFrame") %>%
    writeOGR(
      dsn= shp_grid_file, 
      layer=".", 
      driver="ESRI Shapefile",
      overwrite=T
    )
  
  raster_file_reescaled_countour <- tempfile() %>% paste0(".tif")
  raster_reescaled_countour <- gdalwarp(
    raster_tmp_file,
    raster_file_reescaled_countour,
    s_srs = raster::crs(raster_stack),
    t_srs = raster::crs(shp),
    cutline = shp_countour_file,
    crop_to_cutline = T,
    r = 'average',
    tr = c(cell_width, cell_height),
    tap = T,
    te = shp %>% bbox() %>% as("vector"),
    te_srs = raster::crs(shp),
    dstnodata = "-9999.0",
    ot = 'Float32',
    co = c("BIGTIFF=YES"), #"COMPRESS=DEFLATE", "PREDICTOR=2","ZLEVEL=9"),
    #wo = c("CUTLINE_ALL_TOUCHED=TRUE"),
    multi = T,
    output_Raster = T,
    overwrite = T,
    verbose = T
  ) %>% 
    raster::crop(shp)

  raster_reescaled_countour_masked <- raster_reescaled_countour %>% 
    terra::rast() 
  
  raster_reescaled_countour_masked <- raster_reescaled_countour_masked %>%
    terra::mask(shp %>% terra::vect(), touches=F) %>% 
    raster::stack()

  raster_grid <- gdal_rasterize(
    shp_area_file,
    tempfile() %>% paste0(".tif"),
    #burn = 0,
    a = "cell_id",
    #at = T, 
    co = c("BIGTIFF=YES"), 
    a_nodata = "-9999.0",    
    tr = c(cell_width, cell_height),
    #tap= T,
    ot = 'Float32',
    output_Raster = T,
    te = shp %>% bbox() %>% as("vector"),
    verbose = T
  )
  
  raster_grid <- raster_grid %>% 
    terra::rast() %>%
    terra::mask((raster_reescaled_countour_masked %>% terra::rast())[[1]])


  grid_cells <- raster_grid %>%
    as.vector() %>%
    discard(is.na)

  shp@data <- shp@data %>%
    as.data.frame()
  
  shp_grid <- shp[grid_cells %>% as.integer(),] 

  shp_grid@data$cell_id <- 1:length(grid_cells)
  shp_grid %>% 
    spChFIDs(as.character(1:length(grid_cells)))


  shp_grid@data <- shp_grid@data %>% bind_cols(
    raster_reescaled_countour_masked %>%
      as.list() %>%
      map_dfc(~ .x %>% values() %>% discard(is.na) %>% as.data.frame()) %>%
      rename_all(~ var_names)
  )
  
  return(shp_grid)
}








