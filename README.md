# scriptsdm





# Preparação do ambiente para execução do script

Obs: Se vc estiver usando as workstations do nupelia, deve primeiro abrir o rstudio, digitar no console o comando: options(repos = c(CRAN = "https://mran.revolutionanalytics.com/snapshot/2020-04-24")); e, atualizar todos os pacotes. Se o RStudio pedir a qualquer momento para reiniciar, não reinicie.


1. Baixar o maxent.jar e colocá-lo na mesma pasta deste documento
   * Baixar de (http://biodiversityinformatics.amnh.org/open_source/maxent/./maxent.php?op=download)
   * Depois da primeira vez que esse script for executado em sua totalidade, esse arquivo pode ser apagado, se você desejar
2. Baixar e instalar o R
3. Baixar e instalar o RTools 3.5
   * Baixar de  (https://cran.r-project.org/bin/windows/Rtools/Rtools35.exe)
   * Ao instalar, marcar a opção para colocar o RTools no *path*
4. Baixar e instalar o RStudio
   * Baixar de (https://download1.rstudio.org/desktop/windows/RStudio-1.2.5001.exe)
5. Instalar e carregar o pacote devtools (execute os comandos usando o console do RStudio)
   * `install.packages("devtools")`
   * `library(devtools)`
6. Instalar o pacote sdm e os pacotes dos quais ele depende (execute os comandos usando o console do RStudio)
   * `install.packages("sdm")`
   * `library(sdm)`
   * `installAll()`
7. Instalar todos os outros pacotes necessários (execute os comandos usando o console do RStudio)
   * `install.packages(c("mclust", "ClusterR", "isotree", "rJava", "lubridate", "snakecase", "janitor", "vroom", "fs", "here", "purrrlyr", "stringr", "dplyr", "purrr", "tidyr", "tibble", "clusternor", "usdm", "rdist", "FactoMineR",   "sdm", "snow", "DT", "Rtsne", "paran", "factoextra", "scales", "ggfortify", "mapview", "plotly",         "cowplot", "patchwork", "ggplot2", "gdalUtils", "rgeos", "sf", "rasterDT", "data.table", "rgdal", "raster", "sp", "parallel", "prettyGraphs", "stars")`          

