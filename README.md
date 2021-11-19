# scriptsdm





# Preparação do ambiente para execução do script

Obs: Se vc estiver usando as workstations **Windows** do nupelia, deve primeiro abrir o rstudio, digitar no console o comando: options(repos = c(CRAN = "https://mran.revolutionanalytics.com/snapshot/2020-04-24")); e, atualizar todos os pacotes. Se o RStudio pedir a qualquer momento para reiniciar, não reinicie.


1. Baixar o maxent.jar e colocá-lo na mesma pasta deste documento
   * Baixar de (http://biodiversityinformatics.amnh.org/open_source/maxent/./maxent.php?op=download)
   * Depois da primeira vez que esse script for executado em sua totalidade, esse arquivo pode ser apagado, se você desejar
2. Baixar e instalar o R
3. Se estiver usando **Windows**, nas workstations do nupélia ou computador pessoal, baixar e instalar o RTools 3.5, caso contrário, ignore esse passo
   * Baixar de  (https://cran.r-project.org/bin/windows/Rtools/Rtools35.exe)
   * Ao instalar, marcar a opção para colocar o RTools no *path*
4. Baixar e instalar o RStudio
   * Baixar de (https://download1.rstudio.org/desktop/windows/)
5. Instalar e carregar o pacote devtools (execute os comandos usando o console do RStudio)
   * `install.packages("devtools")`
   * `library(devtools)`
6. Instalar o pacote sdm e os pacotes dos quais ele depende (execute os comandos usando o console do RStudio)
   * `install.packages("sdm")`
   * `library(sdm)`
   * `installAll()`
7. Instalar todos os outros pacotes necessários (execute os comandos usando o console do RStudio)
   * `install.packages(c("raster", "terra", "rgdal", "rasterDT", "sf", "rgeos", "gdalUtils", "stars", "maptools", "magrittr", "ggplot2", "patchwork", "cowplot", "plotly", "mapview", "ggfortify", "scales", "factoextra", "paran", "Rtsne", "ggcorrplot", "DT", "parallel", "snow", "sdm", "FactoMineR", "rdist", "usdm", "ade4", "tidyverse", "purrrlyr", "here", "fs", "stringr", "vroom", "janitor", "snakecase", "lubridate", "stars", "sf", "ncdf4", "tidync")`          

