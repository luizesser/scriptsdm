library(stars)
library(here)
library(tidyverse)

nome_shape_area_estudo <- here("/home/reginaldo/UTFPR/projetos/modelagem_de_especies/castela_tweedii/entrada/shape_area_estudo/lim_pais_a[6933].shp")  


shape_area_estudo <- sf::read_sf(nome_shape_area_estudo, quiet = T)

repair_shp <- function(shp) {
  #if (sf::st_is(x, c("POLYGON", "MULTIPOLYGON")))
  #if (!(sf::st_geometry(shp) %>% inherits(c("sfc_POLYGON","sfc")))){
  if (!(sf::st_geometry(shp) %>% class() %in% "sfc" %>% any())){
    stop("Arquivo da área de estudo inválido!")
  }
  if (shp %>% sf::st_crs() %>% is.na()){
    shp <- shp %>% 
      sf::st_set_crs(4326)
  } 
  shp <- shp %>% sf::st_make_valid() 
  return(shp)
}

drop_areas <- function(shp, inf_limit=0) {
  if (!(sf::st_geometry(shp) %>% class() %in% "sfc" %>% any())){
    stop("Arquivo da área de estudo inválido!")
  }
  if (inf_limit<=0){
    stop("Limite de área de inválido!")
  }
  
  if (shp %>% sf::st_is(c("POLYGON", "MULTIPOLYGON")) %>% all()) {
    kept_areas <- shp[shp %>% sf::st_area() >= units::set_units(inf_limit,"m^2"),]

    if (shp %>% sf::st_area() %>% sum() != kept_areas %>% sf::st_area() %>% sum()) {
      return(kept_areas)
    }
  }
  return(shp)      
}

make_grid <- function(shp, cell_width=0, cell_height=0, var_names=NULL, centroid=T){
  sfc_as_cols <- function(x, geometry, names = c("x","y")) {
    if (missing(geometry)) {
      geometry <- sf::st_geometry(x)
    } else {
      geometry <- rlang::eval_tidy(enquo(geometry), x)
    }
    stopifnot(inherits(x,"sf") && inherits(geometry,"sfc_POINT"))
    ret <- sf::st_coordinates(geometry)
    ret <- tibble::as_tibble(ret)
    stopifnot(length(names) == ncol(ret))
    x <- x[ , !names(x) %in% names]
    ret <- setNames(ret,names)
    dplyr::bind_cols(x,ret)
  }
  
  if (!(sf::st_geometry(shp) %>% class() %in% "sfc" %>% any())){
    stop("Arquivo da área de estudo inválido!")
  }
  if (cell_width<=0 || cell_height<=0){
    stop("Tamanho de célula inválido!")
  }
  shp <- shp %>%
    dplyr::rename_all(tolower)
  
  if (!is.null(var_names)){
    shp <- shp %>%
      dplyr::select(var_names %>% tolower() %>% all_of()) 
  }
  
  bbox <- shp %>% 
    sf::st_bbox() %>% 
    as.vector()
  
  xmin <- bbox[1]
  ymin <- bbox[2]
  xmax <- bbox[3]
  ymax <- bbox[4]
  
  n_x_cell <- (abs(xmax - xmin) / cell_width) %>%
    ceiling()
  n_y_cell <- (abs(ymax - ymin) / cell_height) %>%
    ceiling()
  
  new_xmin <- ((xmin / cell_width) %>% trunc()) * cell_width
  new_ymin <- ((ymin / cell_height) %>% trunc()) * cell_height
  new_xmax <- new_xmin + (n_x_cell * cell_width)
  new_ymax <- new_ymin + (n_y_cell * cell_height)
  
  attr(sf::st_geometry(shp), "bbox") <- sf::st_bbox(
    c(xmin=new_xmin, xmax=new_xmax, ymax=new_ymax, ymin=new_ymin) %>% 
      sf::st_bbox()
  )
  
  shp_grid <- shp %>% 
    sf::st_make_grid(what = "polygons", n=c(n_x_cell, n_y_cell))
  
  shp_grid <- st_sf(geometry=shp_grid) %>%
    rowid_to_column("cell_id")
  
  
  grid_cells <- shp_grid %>% 
    sf::st_intersects(shp) %>% 
    as.data.frame() %>%
    rename(cell_id=row.id, area_id=col.id) %>%
    dplyr::left_join(shp_grid, by = "cell_id")
  
  shp_grid <- grid_cells %>% 
    dplyr::left_join(
      shp %>% 
        as.data.frame() %>% 
        tibble::rowid_to_column("area_id") %>% 
        dplyr::select(-geometry),
      by = "area_id"
    ) %>%
    dplyr::group_by(cell_id) %>% 
    dplyr::summarise_all(~ ifelse(is.numeric(.), mean(.), .), na.rm = TRUE) %>%
    dplyr::select(-area_id, -cell_id) %>%
    tibble::rowid_to_column("cell_id") %>%
    st_sf(crs=sf::st_crs(shp), sf_column_name = "geometry")
  
  if (centroid){
    centroids <-  shp_grid %>%
      sf::st_geometry() %>%
      sf::st_centroid() %>% 
      sf::st_coordinates() %>%
      as.data.frame()

    shp_grid <- shp_grid %>%
      dplyr::bind_cols(centroids) %>%
      dplyr::rename(x_centroid = X, y_centroid = Y)
  }
  return(shp_grid)
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


