FROM openanalytics/r-base:4.0.5
MAINTAINER Luis V. Valcarcel "lvalcarcelg@gmail.com" 

# System libraries of general use

RUN apt-get update && apt-get upgrade -y binutils && apt-get clean 

RUN apt-get update \
    && apt-get install -y \
    sudo \
    pandoc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libv8-dev \
    libgsl-dev \
    libssh2-1-dev \
    curl \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean 


RUN apt-get update --allow-insecure-repositories --allow-unauthenticated && apt-get install -y --allow-unauthenticated  \
  libxml2-dev \
  libudunits2-dev \
  libssl-dev \
  libfontconfig1-dev \
  libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev 
  # && apt-get clean && apt-get autoremove -y
  #  libcurl4-gnutls-dev \
#  libcurl4-openssl-dev \
# r-base-core

RUN apt update && apt install -y cmake
		
COPY Rprofile.site /usr/lib/R/etc/

# Install dependencies of the app
#RUN R -e "pkgsNeeded <- c('rgeos', 'rgdal'); install.packages(pkgsNeeded, dependencies = T, repos = 'https://cloud.r-project.org'); q()" 
#RUN R -e "pkgsNeeded <- c('htmltools', 'nloptr'); install.packages(pkgsNeeded, dependencies = T, repos = 'https://cloud.r-project.org'); q()" 

RUN R -e "pkgsNeeded <- 'https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.5.1.tar.gz'; install.packages(pkgs = pkgsNeeded, repos = NULL); q()" 
RUN R -e "pkgsNeeded <- c('shiny','shinyjs','shinythemes','shinyWidgets', 'rmarkdown','shinycssloaders', 'shinydashboard', 'shinydashboardPlus','shinyalert', 'shinyhelper', 'DT', 'shinyTree'); pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages()); if (length(pkgsNeeded[!pkgsInstalled]) > 0){print(paste(pkgsNeeded[!pkgsInstalled])); install.packages(pkgsNeeded[!pkgsInstalled], dependencies = T, repos = 'https://cloud.r-project.org')}; q()" 
RUN R -e "pkgsNeeded <- c('tidyverse'); pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages()); if (length(pkgsNeeded[!pkgsInstalled]) > 0){print(paste(pkgsNeeded[!pkgsInstalled])); install.packages(pkgsNeeded[!pkgsInstalled], dependencies = T, repos = 'https://cloud.r-project.org')}; q()" 

RUN R -e "pkgsNeeded <- c('data.table', 'dplyr', 'rstatix','Matrix', 'reshape2', 'gridExtra', 'pheatmap','ggplot2', 'ggpubr', 'ggrepel', 'scales', 'colorspace', 'circlize', 'BiocManager', 'openxlsx', 'viridis', 'latex2exp', 'benchmarkme'); pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages()); if (length(pkgsNeeded[!pkgsInstalled]) > 0){print(paste(pkgsNeeded[!pkgsInstalled])); install.packages(pkgsNeeded[!pkgsInstalled], dependencies = T, repos = 'https://cloud.r-project.org')}; q()" 
RUN R -e "pkgsNeeded <- c('foreach', 'doParallel', 'parallel', 'doSNOW', 'doParallel', 'snow', 'pbapply'); pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages()); if (length(pkgsNeeded[!pkgsInstalled]) > 0){print(paste(pkgsNeeded[!pkgsInstalled])); install.packages(pkgsNeeded[!pkgsInstalled], dependencies = T, repos = 'https://cloud.r-project.org')}; q()" 
RUN R -e "pkgsNeeded <- c('Rgraphviz', 'RBGL'); pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages()); if (length(pkgsNeeded[!pkgsInstalled]) > 0){print(paste(pkgsNeeded[!pkgsInstalled])); BiocManager::install(pkgsNeeded[!pkgsInstalled], ask = F)}; q()" 
RUN R -e "pkgsNeeded <- c('RcppZiggurat', 'Rfast'); pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages()); if (length(pkgsNeeded[!pkgsInstalled]) > 0){print(paste(pkgsNeeded[!pkgsInstalled])); install.packages(pkgsNeeded[!pkgsInstalled], dependencies = T, repos = 'https://cloud.r-project.org')}; q()" 
RUN R -e "pkgsNeeded <- c('gRbase', 'Rfast2', 'parallel', 'ParallelLogger','pbmcapply', 'doParallel', 'foreach', 'parallelly', 'future', 'future.apply'); pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages()); if (length(pkgsNeeded[!pkgsInstalled]) > 0){print(paste(pkgsNeeded[!pkgsInstalled])); install.packages(pkgsNeeded[!pkgsInstalled], dependencies = T, repos = 'https://cloud.r-project.org')}; q()" 


RUN R -e "pkgsNeeded <- c('circlize', 'colorspace', 'data.table', 'doParallel', 'doSNOW', 'dplyr', 'DT', 'foreach', 'ggplot2', 'ggpubr', 'ggrepel', 'gridExtra', 'latex2exp', 'Matrix', 'openxlsx', 'parallel', 'pbapply', 'pheatmap', 'reshape2', 'rmarkdown', 'rstatix', 'scales', 'shiny', 'shinyalert', 'shinycssloaders', 'shinydashboard', 'shinydashboardPlus', 'shinyhelper', 'shinyjs', 'shinythemes', 'shinyTree', 'shinyWidgets', 'snow', 'tidyverse', 'viridis'); pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages()); if (length(pkgsNeeded[!pkgsInstalled]) > 0){print(paste(pkgsNeeded[!pkgsInstalled])); install.packages(pkgsNeeded[!pkgsInstalled], dependencies = T, repos = 'https://cloud.r-project.org')}; q()" 

#RUN R -e "pkgsNeeded <- c('abind', 'backports', 'base64enc', 'boot', 'brio', 'broom', 'bslib', 'cachem', 'callr', 'car', 'carData', 'circlize', 'cli', 'codetools', 'colorspace', 'commonmark', 'corrplot', 'cowplot', 'cpp11', 'crayon', 'crosstalk', 'data.table', 'desc', 'diffobj', 'digest', 'doParallel', 'doSNOW', 'dplyr', 'DT', 'ellipsis', 'evaluate', 'fansi', 'farver', 'fastmap', 'fontawesome', 'foreach', 'foreign', 'fresh', 'fs', 'generics', 'ggplot2', 'ggpubr', 'ggrepel', 'ggsci', 'ggsignif', 'GlobalOptions', 'glue', 'graphics', 'grDevices', 'grid', 'gridExtra', 'gtable', 'highr', 'htmltools', 'htmlwidgets', 'httpuv', 'isoband', 'iterators', 'jquerylib', 'jsonlite', 'knitr', 'labeling', 'later', 'latex2exp', 'lattice', 'lazyeval', 'lifecycle', 'lme4', 'magrittr', 'maptools', 'markdown', 'MASS', 'Matrix', 'MatrixModels', 'methods', 'mgcv', 'mime', 'minqa', 'munsell', 'nlme', 'nloptr', 'nnet', 'numDeriv', 'openxlsx', 'parallel', 'pbapply', 'pbkrtest', 'pheatmap', 'pillar', 'pkgconfig', 'pkgload', 'plyr', 'polynom', 'praise', 'processx', 'promises', 'ps', 'purrr', 'quantreg', 'R6', 'rappdirs', 'RColorBrewer', 'Rcpp', 'RcppEigen', 'rematch2', 'reshape2', 'rlang', 'rmarkdown', 'rprojroot', 'rstatix', 'rstudioapi', 'sass', 'scales', 'shape', 'shiny', 'shinyalert', 'shinycssloaders', 'shinydashboard', 'shinydashboardPlus', 'shinyhelper', 'shinyjs', 'shinythemes', 'shinyTree', 'shinyWidgets', 'snow', 'sourcetools', 'sp', 'SparseM', 'splines', 'stats', 'stringi', 'stringr', 'testthat', 'tibble', 'tidyr', 'tidyselect', 'tinytex', 'tools', 'tydyverse', 'utf8', 'utils', 'uuid', 'vctrs', 'viridis', 'viridisLite', 'waiter', 'waldo', 'withr', 'xfun', 'xtable', 'yaml', 'zip'); pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages()); if (length(pkgsNeeded[!pkgsInstalled]) > 0){print(paste(pkgsNeeded[!pkgsInstalled])); install.packages(pkgsNeeded[!pkgsInstalled], dependencies = T, repos = 'https://cloud.r-project.org')}; q()" 
#RUN R -e "pkgsNeeded <- c('limma', 'tximport'); pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages()); if (length(pkgsNeeded[!pkgsInstalled]) > 0){print(paste(pkgsNeeded[!pkgsInstalled])); BiocManager::install(pkgsNeeded[!pkgsInstalled])};q()" 
# RUN R -e "pkgsNeeded <- c('benchmarkme', 'shinyTree','shinyjs','shinythemes','shinyWidgets', 'rmarkdown', 'shinycssloaders', 'shinydashboardPlus','shinyalert', 'shinyhelper', 'DT', 'Matrix', 'foreach', 'doParallel', 'parallel', 'data.table', 'doSNOW', 'reshape2', 'ggplot2', 'gridExtra', 'pheatmap', 'ggpubr', 'openxlsx', 'BiocManager', 'viridis', 'latex2exp', 'curl', 'rio', 'car', 'rstatix', 'pbapply', 'feather'); pkgsInstalled <-  pkgsNeeded %in% rownames(installed.packages()); if (length(pkgsNeeded[!pkgsInstalled]) > 0){print(paste(pkgsNeeded[!pkgsInstalled])); install.packages(pkgsNeeded[!pkgsInstalled], dependencies = T, repos = 'https://cloud.r-project.org')}; q()" 

EXPOSE 3838

RUN apt-get update && apt-get install -y git tree

# Copy the app to the image
RUN mkdir /root/gmctool 
#RUN tree /root/

# trick to skip cache
#ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
#RUN git clone --depth 1 --progress https://github.com/lvalcarcel/gmctool.git /root/github && cp -rv /root/github/shinyApp/* /root/gmctool &&  rm -rvf /root/github

ADD "https://api.github.com/repos/lvalcarcel/gmctool/commits?per_page=1" latest_commit
RUN curl -L "https://github.com/lvalcarcel/gmctool/archive/main.zip" -o /root/main.zip && cd /root/ && unzip /root/main.zip && cp -rv /root/gmctool-main/shinyApp/* /root/gmctool &&  rm -rvf /root/gmctool-main/ && rm -v  /root/main.zip

#COPY gMCStool /root/gMCStool

# comment for more powerfull cluster
COPY configuration_computer_power.R /root/gmctool/configuration_computer_power.R 
RUN tree /root/


RUN chmod 777 /root/gmctool/www/* && chmod 777 /root/gmctool/instructions/* && cd /root/gmctool/

# RUN cd /root/gMCStool/
# CMD ["R", "-e", "shiny::runApp('/root/gMCStool')"]

