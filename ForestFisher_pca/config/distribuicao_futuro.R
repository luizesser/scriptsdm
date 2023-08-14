summ_sp_predictions_from_folder <- function(pred_folder){
  folders <- pred_folder %>% path("predicao_sp") %>%
    dir_ls(type = "directory")
  
  if (pred_folder %>% path("consenso_sp") %>% dir_exists() == F){
    pred_folder %>% 
      path("consenso_sp") %>%
      dir_create()
  }
  
  for (folder in folders) {
    sp <- folder %>% 
      path_file()
    if (pred_folder %>% path("consenso_sp") %>% path(sp) %>% dir_exists() == F){
      pred_folder %>% 
        path("consenso_sp") %>% 
        path(sp) %>%
        dir_create()
      
      suppressMessages(
        df_pred_freq <- folder %>% 
          dir_ls(type = "file", recurse = T, glob = "*_freq.csv") %>% 
          purrr::map(~ vroom(.x) %>% select(consenso)) %>% 
          purrr::set_names(~ str_remove(path_file(.x), "_freq.csv")) %>% 
          purrr::imap(~ set_names(.x, .y)) %>% 
          map_dfc(~ (.)) 
      )
      
      df_pred_freq <-rep(sp, nrow(df_pred_freq)) %>%
        as.data.frame() %>%
        rename("especie"=".") %>%
        bind_cols(df_pred_freq)
      
      suppressMessages(
        df_pred_pa <- folder %>% 
          dir_ls(type = "file", recurse = T, glob = "*_pa.csv") %>% 
          purrr::map(~ vroom(.x) %>% select(consenso)) %>% 
          purrr::set_names(~ str_remove(path_file(.), "_pa.csv")) %>% 
          purrr::imap(~ set_names(.x, .y)) %>% 
          map_dfc(~ (.))
      )
      
      df_pred_pa <- rep(sp, nrow(df_pred_pa)) %>%
        as.data.frame() %>%
        rename("especie"=".") %>%
        bind_cols(df_pred_pa)
      
      suppressMessages(
        df_pred_freq %>%
          vroom_write(
            pred_folder %>% path("consenso_sp") %>% path(sp) %>% path(paste0(sp, "_freq.csv")),
            delim = ";"
          )
      )
      suppressMessages(
        df_pred_pa %>%
          vroom_write(
            pred_folder %>% path("consenso_sp") %>% path(sp) %>% path(paste0(sp, "_pa.csv")),
            delim = ";"
          )
      )
    }
  }
}

summ_scenarios_predictions_from_folder <- function(pred_folder){
  if (pred_folder %>% path("consenso_cenarios") %>% dir_exists() == F){
    pred_folder %>% 
      path("consenso_cenarios") %>% 
      dir_create()
  }
  
  suppressMessages(
    df_pred <- pred_folder %>% 
        path("consenso_sp") %>% 
        dir_ls(type = "file", glob = "*_freq.csv", recurse = T)  %>%
        vroom()
  )
  
  sp_names <- df_pred$especie %>% 
    unique()
  
  scenarios_names <- df_pred %>% 
    colnames() %>% 
    purrr::discard(. == "especie")
  
  number_of_rows <- df_pred %>% nrow() / sp_names %>% length()
  
  for (scenario in scenarios_names){
    if (pred_folder %>% path("consenso_cenarios") %>% path(scenario) %>% dir_exists() == F){
      pred_folder %>% 
        path("consenso_cenarios") %>% 
        path(scenario) %>% 
        dir_create()
      
      suppressMessages(
        df_pred[[scenario]] %>% 
          matrix(nrow=number_of_rows) %>% 
          as.data.frame() %>% 
          set_names(sp_names) %>%
          vroom_write(
            pred_folder %>% path("consenso_cenarios") %>% path(scenario) %>% path(paste0(scenario, "_freq.csv")),
            delim = ";"
          )
      )
    }  
  }
  
  suppressMessages(
    df_pred <- pred_folder %>% 
      path("consenso_sp") %>% 
      dir_ls(type = "file", glob = "*_pa.csv", recurse = T)  %>%
      vroom()
  )
  
  for (scenario in scenarios_names){
    suppressMessages(
      df_pred[[scenario]] %>% 
        matrix(nrow=number_of_rows) %>% 
        as.data.frame() %>% 
        set_names(sp_names) %>%
        vroom_write(
          pred_folder %>% path("consenso_cenarios") %>% path(scenario) %>% path(paste0(scenario, "_pa.csv")),
          delim = ";"
        )
    )
  }
}

