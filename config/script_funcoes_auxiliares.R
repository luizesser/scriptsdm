# Aqui são carregadas as bibliotecas e configurações necessárias para a execução dos documentos
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE, fig.width = 12, out.width="100%")

# Geoprocessamento
library(raster)                              # manipulação de rasters
#library(terra)                              # manipulação de rasters
#detach("package:terra", unload=T)
library(rgdal)                               # readOGR, ler arquivos raster e shapefile
library(rasterDT)                            # rasterização
library(sf)                                  # manipulação de formas geométricas com interface com gdal, rgeos e proj4g
library(rgeos)                               # habilita operações geométricas
library(gdalUtils)                           #  
library(stars)                               #
library(maptools)
library(magrittr)

# Gráficos
library(ggplot2)                             # plotagem de gráficos
library(patchwork)                           # combinar gráficos ggplot independentes
library(cowplot)                             # combinar gráficos ggplot independentes
library(plotly)                              # adiciona interatividade nos gráficos
library(mapview)                             # ferramenta de visualização de mapas geográficos
library(ggfortify)                           # permite que o ggplot plote gráficos rasters e shapefiles
library(scales)                              # habilita o uso de diferentes projeções gegráficas nos gráficos de mapas
library(factoextra)                          # Visualização de PCA
library(paran)                               # Análise paralela de Horn
library(Rtsne)                               # Redutor de dimensionalidade T-SNE  
library(ggcorrplot)

# Tabelas
library(DT)                                  # renderiza tabelas interativas

# Processamento e computação paralela 
library(parallel)                            # habilita processamento paralelo
library(snow)                                # habilita computação paralela

# SDM
library(sdm)                                 # modelagem de distribuição de espécies
library(FactoMineR)                          # PCA
library(rdist)                               # distancia euclideana
library(usdm)                                # vif
#library(clusternor)                          # xmeans clustering
library(ade4)                                # Enfa
library(boot)                                # Inverse Logit Function for glmnet

# Miscelânea
library(tidyverse)                           # funções de manipulação de dataframes e  listas
library(purrrlyr)                            # funções que auxiliam no uso conjunto do dplyr com o purrr
library(here)                                # encontra o caminho (pasta) do projeto de maneira inteligente  e segura
library(fs)                                  # funções de manipulação de pastas e arquivos
library(stringr)                             # manipulação de strings
library(vroom)                               # le arquivos csv
library(janitor)                             # faz limpeza de datasets e nomes de variáveis
library(snakecase)                           # transformação de strings em snake case
library(lubridate)                           # manipulação de datas
library(rgbif)                               # download GBIF data.

select <- dplyr::select
here <- here::here
map <- purrr::map

options(java.parameters = "-Xmx1g", java.awt.headless="true")
source(here("config/ggcorrplot.R"))
source(here("config/xmeans.R"))
source(here("config/algoritmos_predicao.R"))
source(here("config/geoprocessamento_novo.R"))
source(here("config/ocorrencias.R"))
source(here("config/avaliacao_variaveis_preditoras.R"))
source(here("config/treinamento_avaliacao.R"))
source(here("config/distribuicao_futuro.R"))


pastaRaiz <- here()

if (interactive()){
  #pastaRaiz <- params$title %>% make_clean_names() %>% here()
} else {
  knitr::opts_chunk$set(echo=TRUE, message=TRUE, warning=TRUE, error=TRUE, include=TRUE)
  myKnit = function(inputFile, encoding, out=NULL){
    #   if (is.null(out)){
    #     out <- inputFile %>% fs::path_file() %>% fs::path_ext_remove()
    #   } else{
    #     out <- out %>% fs::path(inputFile %>% fs::path_file() %>% fs::path_ext_remove())
    #   }
    # pastaRaiz <<- out
    # rmarkdown::render(inputFile, encoding = encoding, output_dir = out)
    rmarkdown::render(inputFile, encoding = encoding, envir = new.env())
  }
  pastaRaiz <- ""
}


#---------------------------------------------------------------------------------------------
# Executa a predição a partir dos modelos treinados e transforma dados de suitability em dados de frequencia e, posteriormente, em presença/ausência

# Essa função recebe como entrada: um dataframe em que linhas são as células da grid e as colunas as varáveis 
# preditoras; os modelos treinados; e uma lista com os nomes dos algoritmos de predição, disponíveis dentro dos modelos treinados,  para fazer a predição.
# A predição é feita para cada célula da grid com cada um dos diversos modelos dos diversos algoritmos de predição desejados. 
# Um limiar é usado para determinar se, em uma dada célula, existe presença ou ausência. O resultado é que, haverá uma matriz
# M x N, onde M são as células e N são as predições dos diversos modelos dos diversos algoritmos. Cada célula terá valor 0 ou 1,
# indicando presença ou ausência na célula. Como podem (e muito provavelmente) haverão diversos modelos de um mesmo algoritmo 
# (em função do uso de replicações como cross-validation ou bootstraping) que indicará ausência para algumas execuções e presença
# para outras execuções, é necessário fazer um consenso para cada um dos algoritmos em cada célula. O consenso é a média dos 
# valores da célula para cada um dos algoritmos. Assim, o resultado é uma matriz M x N, onde M são as células da grid e N são os 
# consensos (médias para aquele algoritmo). O valor de cada célula será um valor entre 0 e 1. 

#df_p = scenarios_list %>% 
#  pluck(scenario_name)
#m_treinados = t_models %>% pluck(sp)
#algoritmo_predicao = pred_methods
#tipo_thresh=2
#lista_thresh=NULL

DRE_predict <- function(df_p, m_treinados, algoritmo_predicao, tipo_thresh=2, lista_thresh=NULL){
  if (is.null(lista_thresh)){
    lista_thresh <- getEvaluation(m_treinados, stat=c('threshold'), wtest="dep.test", opt=tipo_thresh)
  }
  predi <- predict(m_treinados, df_p, method=algoritmo_predicao) %>% 
    as.data.frame() %>%
    replace(., is.na(.), 0)
  
  predi[is.na(predi[])] <- 0
  
  colnames(predi) <- colnames(predi) %>% 
    modify(., ~ str_remove(unlist(str_split(., "[.-]"))[1], "id_"))
  
  for (nome_coluna in colnames(predi)){
    thresh <- lista_thresh %>% 
      filter(modelID==as.integer(nome_coluna))
    
    thresh <- thresh$threshold
    if (length(thresh)>0){
      predi[predi[nome_coluna] >= thresh, nome_coluna] <- 1
      predi[predi[nome_coluna] < thresh, nome_coluna] <- 0
    }
  }
  
  detalhes_modelos <- getModelInfo(m_treinados) %>% 
    filter(method %in% algoritmo_predicao) %>%
    select(modelID, method)
  
  for (i in 1:ncol(predi)){
    colnames(predi)[i] <- as.character(detalhes_modelos$method[detalhes_modelos$modelID==as.integer(colnames(predi[i]))])
  }
  
  colunas_por_metodo <- ncol(predi) / length(unique(colnames(predi)))
  predi_temp = data.frame(row.names = 1:nrow(predi))
  for (i in 1:length(unique(colnames(predi)))){
    nome_coluna <- colnames(predi)[(i-1)*colunas_por_metodo+1]
    predi_temp[[nome_coluna]] <- predi[,((i-1)*colunas_por_metodo+1):(i*colunas_por_metodo)] %>% rowSums(.)
    predi_temp[[nome_coluna]] <- predi_temp[[nome_coluna]] / colunas_por_metodo
  }
  
  return(predi_temp)
}

# ------------------------------------------------------------------------------------------------------------------------
presences_number <- function(grid_matrix_pa, spp_names){
  spp_names_abv <- spp_names %>% 
    to_snake_case() %>% 
    abbreviate(minlength = 10) %>% 
    as.vector()
  x <- sapply(spp_names_abv, 
              function(x){
                sum(as.data.frame(grid_matrix_pa)[,x])},
              USE.NAMES = T) %>% data.frame(Presences=.)
  return(x)
}

richness_map <- function(df_pred, shp_estudo){
  shp_estudo <- cbind(shp_estudo, df_pred)
  df <- as.data.frame(shp_estudo)
  df$df_pred <- as.numeric(df$df_pred)
  
  mapa_temp <- ggplot(st_as_sf(df)) +
    geom_sf(aes(fill = df_pred), color=NA) +
    scale_fill_continuous(type='viridis', limits=c(0,max(df$df_pred))) +
    ggtitle(paste0('Richness'))
  
  return(mapa_temp)
}

ensemble_map <- function(df_pred, shp_estudo){
  shp_estudo <- cbind(shp_estudo, df_pred)
  df <- as.data.frame(shp_estudo)
  df$df_pred <- as.numeric(df$df_pred)
  
  mapa_temp <- ggplot(st_as_sf(df)) +
      geom_sf(aes(fill = df_pred), color=NA) +
      scale_fill_continuous(low="#D1392C", high="#4A7CB3", limits=c(0,max(df$df_pred))) +
      ggtitle(paste0('Ensemble'))

  return(mapa_temp)
}

distribution_map <- function(df_pred, shp_estudo, returnRasters=F){
  shp_estudo <- cbind(shp_estudo, df_pred)
  df <- as.data.frame(shp_estudo)
  
  if(length(unique(df$consensus)) <= 2){
    df$consensus <- as.factor(df$consensus)
    mapa_temp <- ggplot(st_as_sf(df)) +
      geom_sf(aes(fill = consensus), color=NA) +
      scale_fill_brewer(palette = "Set1")+
      ggtitle(unique(df_pred$species))
  } else {
    mapa_temp <- ggplot(st_as_sf(df)) +
      geom_sf(aes(fill = consensus), color=NA) +
      scale_fill_continuous(low="#D1392C", high="#4A7CB3", limits=c(0,1)) +
      ggtitle(unique(df_pred$species))
  }
  if(returnRasters==T){
    st_rasters <- st_rasterize(shp_estudo %>% dplyr::select("consensus", geometry))
    return (st_rasters)
  } else {
    return(mapa_temp)
  }
}

distribution_map_algo <- function(df_pred, shp_estudo, nome_metodos, returnRasters=F){
  shp_estudo <- cbind(shp_estudo, df_pred)
  nome_metodos <- nome_metodos %>% 
    to_snake_case() %>% 
    abbreviate(minlength = 10) %>%
    unname()
  df <- as.data.frame(shp_estudo)
  if(length(unique(df[,nome_metodos[1]])) <= 2){
    cols <- colnames(df) %in% nome_metodos
    df[,cols] %<>% lapply(function(x){as.factor(x)}) 
    mapa_temp <- lapply(nome_metodos, function(algo){
      ggplot(st_as_sf(df)) +
        geom_sf(aes(fill = !!sym(algo)), color=NA) +
        scale_fill_brewer(palette = "Set1")})
  } else {
    mapa_temp <- lapply(nome_metodos, function(algo){
      ggplot(st_as_sf(df)) +
        geom_sf(aes(fill = !!sym(algo)), color=NA) +
        scale_fill_continuous(low="#D1392C", high="#4A7CB3", limits=c(0,1))})
  }
  if(returnRasters==T){
    st_rasters <- sapply(nome_metodos, function(algo){st_rasterize(shp_estudo %>% dplyr::select(algo, geometry))}, simplify=F, USE.NAMES=TRUE)
    return (st_rasters)
  } else {
    return(mapa_temp)
  }
}

distribution_map_ensembles <- function(df_pred, shp_estudo, ensemble_method=c('wmean_AUC', 'wmean_TSS'), returnRasters=F){ ### Output=lista fazer outra função pra plotar uma imagem com os ids da lista.
  shp_estudo <- cbind(shp_estudo, df_pred)
  ensemble_method <- ensemble_method %>% 
    to_snake_case() %>% 
    abbreviate(minlength = 10) %>%
    unname()
  df <- as.data.frame(shp_estudo)
  
    mapa_temp <- lapply(ensemble_method, function(ens){
      ggplot(st_as_sf(df)) +
        geom_sf(aes(fill = !!sym(ens)), color=NA) +
        scale_fill_continuous(low="#D1392C", high="#4A7CB3", limits=c(0,1))})

  if(returnRasters==T){
    st_rasters <- sapply(nome_metodos, function(algo){st_rasterize(shp_estudo %>% dplyr::select(algo, geometry))}, simplify=F, USE.NAMES=TRUE)
    return (st_rasters)
  } else {
    return(mapa_temp)
  }
}

mapaDistrModelos <- function(df_distr, grid_estudo, nome_modelos){
  return (distribution_map_algo(df_distr, grid_estudo, nome_modelos))
}

matriz_confusao <- function(obs, pre) {
  cmx<-matrix(nrow=2,ncol=2)
  colnames(cmx) <- rownames(cmx) <- c('P','A')
  cmx[1,1] <- length(which(obs == 1 & pre == 1))
  cmx[2,2] <- length(which(obs == 0 & pre == 0))
  cmx[1,2] <- length(which(obs == 0 & pre == 1))
  cmx[2,1] <- length(which(obs == 1 & pre == 0))
  cmx[] <- as.numeric(cmx)
  return (cmx)
}


calcular_metricas_matriz <- function(mc) {
  TP<-mc[1,1]
  FP<-mc[1,2]
  TN<-mc[2,2]
  FN<-mc[2,1]
  valor_mcc <- ((TP*TN)-(FP-FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  valor_recall <- (TP / (TP+FN))
  
  return (round(c(tp=TP, fp=FP, tn=TN, fn=FN, mcc=valor_mcc, recall=valor_recall), 3))
}


calcular_limiares_modelo <- function(obs, pre){
  th <- sort(unique(pre))
  aval <- matrix(nrow=length(th),ncol=7)
  colnames(aval) <- c("pre", "tp","fp","tn","fn","mcc", "recall")
  aval[,1] <- sort(unique(pre))
  
  for (i in seq_along(th)) {
    w <- which(pre >= th[i])
    pt <- rep(0,length(pre))
    pt[w] <- 1
    aval[i,2:7] <- calcular_metricas_matriz(matriz_confusao(obs,pt))
  }
  return (as.data.frame(aval))
}

calcular_limiar_metrica<- function(limiares, tipo="mcc"){
  return ((limiares %>%
             filter(get(tipo) == max(get(tipo))))$pre %>%
            max())
}


calcular_limiar_todos_modelos <- function(modelos){
  aval_list <- list()
  modelos <- modelos@models$especie
  for (alg in modelos){
    for (modelo in alg){
      if (length(modelo@evaluation)<=0){
        aval_list[modelo@mID] <- 0
      }
      else {
        aval_list[modelo@mID] <- calcular_limiares_modelo(modelo@evaluation$test.dep@observed,
                                                          modelo@evaluation$test.dep@predicted) %>%
          calcular_limiar_metrica("mcc")
      }
    }
  }
  
  return (aval_list %>% map_df(~ as.data.frame(.), .id="modelID") %>% rename(threshold_MCC="."))
}


knit_chunk_as_child <- function(){
  if (is.null(knitr::all_labels())){
    cat("You must knit this document!", sep="\n")
    return()
  }
  if (knitr::opts_current$get("label") %>% str_detect("chunk") %>% any()){
    cat("This chunk must not be include as a child document!", sep="\n")
    return()
  }
  
  current_folder <- knitr::opts_current$get("label") %>%
    discard(~ (.) %>% str_detect("chunk") %>% any()) %>% 
    keep(~ (.) %>% here() %>% is_dir()) %>% 
    here() %>% 
    path() %>% 
    dir_info(type = "directory") %>%
    filter(path %>% str_detect("output_data")) %>%
    select(path, modification_time)
  
  if (nrow(current_folder)==0){
    return(
      here(knitr::opts_current$get("label"),paste0(knitr::opts_current$get("label"), ".Rmd")) %>%
        knitr::knit_child(quiet = TRUE) %>%
        cat(sep = '\n')
    )
  } 
  
  all_folders_details <- knitr::all_labels() %>% 
    discard(~ (.) %>% str_detect("chunk") %>% any()) %>%
    keep(~ (.) %>% here() %>% is_dir()) %>% 
    here() %>% 
    path() %>% 
    dir_info(type = "directory") %>%
    filter(path %>% str_detect("output_data")) %>%
    select(path, modification_time) %>%
    rownames_to_column(var = "id") %>%
    mutate(time_dif = as_datetime(.$modification_time) - lag(as_datetime(.$modification_time))) %>%
    replace(is.na(.), 0)
  
  
  is_newest <- all_folders_details %>%
     filter(path == current_folder$path) %>%
     select(time_dif) %>%
     pull() %>%
    {(.<0)}
  
  if (is_newest){
    return(
      here(knitr::opts_current$get("label"),paste0(knitr::opts_current$get("label"), ".Rmd")) %>%
        knitr::knit_child(quiet = TRUE) %>%
        cat(sep = '\n')
    )
  }    
  
  is_last_one <- all_folders_details %>%
    filter(path == current_folder$path) %>%
    select(id) %>%
    pull() %>%
    {(. == all_folders_details$id %>% last())}
  
  if (is_last_one)
    all_folders_details %>%
      datatable(options = list(pageLength = 10, scrollX=T))
}


render_chunk_as_html <- function(){
  current_chunk <- knitr::opts_current$get("label")
  all_chunks <- knitr::all_labels()

  if (is.null(all_chunks)){
    cat("You must knit this document!", sep="\n")
    return()
  }
  if (current_chunk %>% str_detect("chunk") %>% any()){
    cat("This chunk must not be include as a child document!", sep="\n")
    return()
  }

  current_folder <- current_chunk %>%
    discard(~ (.) %>% str_detect("chunk") %>% any()) %>% 
    keep(~ (.) %>% here() %>% is_dir()) %>% 
    here() %>% 
    path() %>% 
    dir_info(type = "directory") %>%
    filter(path %>% str_detect("output_data")) %>%
    select(path, modification_time)
  
  if (nrow(current_folder)==0){
    path <- callr::r(
      function(...) rmarkdown::render(...),
      args = list(
        here(current_chunk,paste0(current_chunk, ".Rmd")), 
        quiet=T, 
        envir = new.env()
      )
    )
  }
  current_folder <- current_chunk %>%
    discard(~ (.) %>% str_detect("chunk") %>% any()) %>% 
    keep(~ (.) %>% here() %>% is_dir()) %>% 
    here() %>% 
    path() %>% 
    dir_info(type = "directory") %>%
    filter(path %>% str_detect("output_data")) %>%
    select(path, modification_time)
  
  all_folders_details <- all_chunks %>% 
    discard(~ (.) %>% str_detect("chunk") %>% any()) %>%
    keep(~ (.) %>% here() %>% is_dir()) %>% 
    here() %>% 
    path() %>% 
    dir_info(type = "directory") %>%
    filter(path %>% str_detect("output_data")) %>%
    select(path, modification_time) %>%
    rownames_to_column(var = "id") %>%
    mutate(time_dif = as.numeric(as_datetime(.$modification_time) - lag(as_datetime(.$modification_time)))) %>% 
    replace(is.na(.), 0)
  
  is_newest <- all_folders_details %>%
    filter(path == current_folder$path) %>%
    select(time_dif) %>%
    pull() %>%
    {(.<0)}
  
  if (is_newest){
    current_folder$path %>%
      dir_delete()
    path <- callr::r(
      function(...) rmarkdown::render(...),
      args = list(
        here(current_chunk,paste0(current_chunk, ".Rmd")), 
        quiet=T, 
        envir = new.env()
      )
    )
  }      
  
  html_src <- current_chunk %>%
    discard(~ (.) %>% str_detect("chunk") %>% any()) %>% 
    keep(~ (.) %>% here() %>% is_dir()) %>% 
    here() %>% 
    path() %>% dir_ls(glob = "*.html") %>% 
    path_rel()
 
  is_last_one <- all_folders_details %>%
    filter(path == current_folder$path) %>%
    select(id) %>%
    pull() %>%
    {(. == all_folders_details$id %>% last())}
  
  if (is_last_one){
    aTable <- all_folders_details %>%
      datatable(filter="none", list(dom = 't')) %>%
      formatDate("modification_time", method="toLocaleString")
    return(
      htmltools::tagList(list(
        htmltools::tags$iframe(
          title = current_chunk %>% make_clean_names(case = "title") ,
          src = html_src,
          width="100%",
          height="900"
        ),
        htmltools::tags$div(
          aTable
        )
      ))
    )
  } else {
    return(
      htmltools::tags$iframe(
        title = current_chunk %>% make_clean_names(case = "title") ,
        src = html_src,
        width="100%",
        height="900"
      )
    )
  }
}

fortify_join <- function(shp){
  shp_temp <- shp %>% 
    fortify(region="cell_id") %>%
    mutate(id = as.integer(id)) %>%
    left_join(shp@data, by=c("id"="cell_id"))
  
  return(shp_temp)
}


### Funções luizfesser:
compare_gcms <- function(folder_future_rasters_gcms, grid_study_area, var_names=c('bio_1','bio_12'), gcm_names, scenarios_gcm, k=NULL){
  # Get stacks files
  l <- list.files(folder_future_rasters_gcms, pattern='.tif', rec=T, full.names = T)
  # Import stacks
  s <- lapply(l[grep(scenarios_gcm,l)], function(x){s <- stack(x)
  names(s) <- paste0('bio_',1:19) # Rename rasters
  return(s)})
  # Name list itens
  names(s) <- sort(gcm_names)
  
  # Transform stacks
  s <- sapply(s, function(x){# Subset stacks to keep only var_names
    x <- x[[var_names]]
    # Reproject to match grid_study_area crs.
    if(!as.character(crs(x))==as.character(CRS(crs(grid_study_area)))){
      x <- projectRaster(x, crs=CRS(crs(grid_study_area)))
    }
    # Crop and mask stacks
    x <- mask(crop(x, grid_study_area),grid_study_area)
    # Transform in data.frames
    x <- x %>% as.data.frame()
    return(x)},
    USE.NAMES = T, 
    simplify = F)
  
  # Scale and flatten variables into one column.
  flatten_vars <- sapply(s, function(x){x <- scale(x)
  x <- as.vector(x)}, USE.NAMES=T)
  
  # Calculate the distance matrix
  dist_matrix <- dist(t(flatten_vars))
  
  # Calculate the means
  m <- sapply(s, function(x){colMeans(x, na.rm=T)}, 
              USE.NAMES = T, simplify = F)
  m <- m %>% as.data.frame() %>% t()
  
  # Create a output
  result <- list()
  
  # Set k
  if(is.null(k)){
    k = 3
  }
  
  # Run K-means
  cl <- kmeans(dist_matrix, k, nstart=1000)
  
  # plot
  kmeans_plot <- fviz_cluster(cl, 
                              data = dist_matrix,
                              palette = "Set1",
                              labelsize = 10,
                              ggtheme = theme_minimal(),
                              main = "K-means Clustering Plot",
                              xlim=c(-3,3),
                              ylim=c(-3,3),
                              legend = 'none')
  
  # Run Hierarchical Cluster
  # hclust_plot <- hclust(dist_matrix)
  # Include elbow, silhouette and gap methods
  flatten_subset <- na.omit(flatten_vars)
  flatten_subset <- flatten_subset[sample(nrow(flatten_subset), nrow(flatten_subset)/20),]
  wss <- fviz_nbclust(flatten_subset, FUN = hcut, method = "wss")
  sil <- fviz_nbclust(flatten_subset, FUN = hcut, method = "silhouette")
  #gap <- fviz_gap_stat(flatten_subset, maxSE = list(method = "globalmax"))
  
  # Compute hierarchical clustering and cut into k clusters
  res <- hcut(t(flatten_subset), k = k)
  dend <- fviz_dend(res, 
                    cex = 0.5, 
                    ylim = c(max(res$height)*1.1/5*-1, max(res$height)*1.1),
                    palette="Set1", 
                    main = "Hierarchical Clustering")
  
  # Run Correlation
  cor_matrix <- cor(flatten_vars, use='complete.obs')
  cor_plot <- ggcorrplot(cor_matrix, 
                         method='circle',
                         type='lower',
                         lab=T, 
                         lab_size = 3,
                         hc.order=T, 
                         hc.method = 'ward.D2', 
                         show.legend = F, 
                         title='Pearson Correlation')
  
  # Plot everything together
  statistics_gcms <- plot_grid(kmeans_plot, 
                               cor_plot, 
                               dend, 
                               plot_grid(wss, 
                                         sil, 
                                         ncol=1, 
                                         labels=c("D", "E"), 
                                         label_x = -0.05),
                               labels = c("A", "B", "C"),
                               ncol = 2, 
                               rel_widths = 4)
  
  gcms <- apply(cl$centers, 1, function(x){which.min(x) %>% names()})
  
  
  return(list(suggested_gcms=as.vector(gcms),
              statistics_gcms=statistics_gcms))
}



predictions_means <- function(predictions_sp, scenarios){
  sp_names <- names(predictions_sp)
  l <- unlist(predictions_sp, recursive = F)
  result <- l[[1]]
  df <- l %>% as.data.frame()
  scenarios <- gsub("-",".",scenarios)
  scenarios2 <- c('current',sort(scenarios))
  for(s in scenarios2){
    result <- cbind(result, rowMeans(df[,grep(paste0(s,'_freq.consensus'), colnames(df))], na.rm = T))
    result <- cbind(result, rowMeans(df[,grep(paste0(s,'_pa.consensus'), colnames(df))], na.rm = T))
    result <- cbind(result, rowSums(df[,grep(paste0(s,'_pa.consensus'), colnames(df))], na.rm = T))
  }
  result <- result[,-c(1:ncol(l[[1]]))]
  names(result) <- sort(c(paste0(scenarios2, '_pa_mean'), paste0(scenarios2, '_pa_sums'), paste0(scenarios2, '_freq_mean')))
  return(result)
}


WorldClim_data <- function(period = 'current', variable = 'bioc', year = '2030', gcm = 'mi', ssp = '126', resolution = 10){
  
  res = ifelse(resolution==30,'s','m')
  
  if(period=='current'){
      if(!dir.exists('input_data/WorldClim_data_current')){ dir.create('input_data/WorldClim_data_current') }
      if(length(list.files("input_data/WorldClim_data_current",pattern='.tif$', full.names=T))==0){
        print(paste0('current_',resolution,res))
        download.file(url = paste0('https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_',
                                   resolution,
                                   res,'_bio.zip'),
                      destfile = paste0('input_data/current_', resolution, res,'.zip'),
                      method = 'auto')
        unzip(zipfile = paste0('input_data/current_', resolution, res,'.zip'),
              exdir = paste0('input_data/WorldClim_data_current'))        
      } else {
        print(paste0('The file for current scenario is already downloaded.'))
      }
   }
  
  if(period=='future'){
    if(!dir.exists('input_data/WorldClim_data_future')){ dir.create('input_data/WorldClim_data_future') }
    all_gcm <- c('ac', 'ae', 'bc', 'ca', 'cc', 'ce','cn', 'ch', 'cr', 'ec','ev', 'fi',
                 'gf', 'gg','gh', 'hg', 'in', 'ic', 'ip', 'me', 'mi', 'mp','ml',
                 'mr', 'uk')
    gcm2 <- c('ACCESS-CM2', 'ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','CanESM5-CanOE','CMCC-ESM2',
              'CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1','EC-Earth3-Veg',
              'EC-Earth3-Veg-LR','FIO-ESM-2-0','GFDL-ESM4','GISS-E2-1-G',
              'GISS-E2-1-H','HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0',
              'IPSL-CM6A-LR','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR',
              'MPI-ESM1-2-LR','MRI-ESM2-0','UKESM1-0-LL')
    gcm3 <- gcm2[match(gcm,all_gcm)]
    all_year <- c('2030', '2050', '2070', '2090')
    year2 <- c('2021-2040', '2041-2060', '2061-2080', '2081-2100')
    year3 <- year2[match(year,all_year)]
    for (g in 1:length(gcm)) {
      for (s in 1:length(ssp)) {
        for (y in 1:length(year)) {
          if(!file.exists(paste0('input_data/WorldClim_data_future/',gcm[g], '_ssp', ssp[s],'_', resolution, '_', year[y],'.tif'))){
            print(paste0(gcm[g], '_ssp', ssp[s], '_', resolution, '_', year[y]))
            download.file(url = paste0('https://geodata.ucdavis.edu/cmip6/',resolution,
                                       res,'/',gcm3[g],'/ssp',ssp[s],'/wc2.1_',resolution,
                                       res,'_',variable,'_',gcm3[g],'_ssp',ssp[s],'_',
                                       year3[y],'.tif'),
                          destfile = paste0('input_data/WorldClim_data_future/',gcm[g], '_ssp', ssp[s],
                                            '_', resolution, '_', year[y],'.tif'),
                          method = 'auto')
          } else {
            print(paste0('The file for future scenario (',
                         paste0('input_data/WorldClim_data_future/',gcm[g], '_ssp', ssp[s],'_', resolution,res, '_', year[y],'.tif'),
                         ') is already downloaded.'))
          }
        }
      }
    }
  }
}


GBIF_data <- function(s, file='input_data/spp_data.csv'){
  if(!file_exists(file)){
    ids <- lapply(s, function(x) { name_suggest(q=x, rank = "species")$data$key[1]})
    ids <- unlist(ids)
    data <- lapply(ids, function(x) { y <- occ_data(taxonKey=x)
    if('decimalLatitude' %in% names(y$data)){
      y <- y$data[,c("species", "decimalLongitude","decimalLatitude")]
      return(y)
    }
    } )
    data <- bind_rows(data)
    data <- data.frame(data)
    data <- na.omit(data)
    data$species <- rep(gsub(' ', '_', s),nrow(data))
    write.csv(data, file, row.names=FALSE)
  } else {
    print(paste0('File already exists. Importing from: ',file))
    data <- read.csv(file)
  }
  return(data)
}


data.clean <- function(x, r=NULL){
  x <- subset( x, !is.na("decimalLongitude") | !is.na("decimalLatitude"))
  x <- cc_cap( x, lon = "decimalLongitude", lat = "decimalLatitude", species = "species")
  x <- cc_cen( x, lon = "decimalLongitude", lat = "decimalLatitude", species = "species")
  x <- cc_dupl(x, lon = "decimalLongitude", lat = "decimalLatitude", species = "species")
  x <- cc_equ( x, lon = "decimalLongitude", lat = "decimalLatitude")
  x <- cc_inst(x, lon = "decimalLongitude", lat = "decimalLatitude", species = "species")
  x <- cc_val( x, lon = "decimalLongitude", lat = "decimalLatitude")
  x <- cc_sea( x, lon = "decimalLongitude", lat = "decimalLatitude")
  if(!is.null(r)){
    print('Raster identified, procceding with rasterized filter.')
    r <- raster(r)
    values(r) <- 1:ncell(r)
    x2 <- x
    coordinates(x2) <- 2:3
    cell_id <- extract(r, x2)
    x <- cbind(x, cell_id)
    x <- x[!duplicated(x[,c(1,4)]),-4]
  }
  return(x)
}


w_area <- function(x){
  cell_size <- area(x, na.rm=TRUE, weights=FALSE)
  x <- cell_size*x
  result <- cellStats(x, sum)
  print(paste0(result," km2"))
  return(as.numeric(result))
}

