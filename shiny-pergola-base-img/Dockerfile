
FROM r-base:3.3.2 

MAINTAINER Jose Espinosa-Carrasco <espinosacarrascoj@gmail.com>

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    libxt-dev \
    libssl-dev \
    #libcurl4-openssl-dev \
    libxml2-dev

# Download and install shiny server
RUN wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb

## Bioconductor dependencies of packages installed from CRAN
# RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("Gviz"); biocLite("GenomicRanges"); biocLite("rtracklayer");'
# RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("GenomicRanges"); biocLite("rtracklayer");'
RUN R -e "install.packages(c('shiny', 'rmarkdown', 'ggplot2', 'XML', 'Rcurl','cowplot', 'dplyr', 'survival', 'gridExtra', 'devtools'), repos='http://cran.rstudio.com/')" \
&& Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("GenomicRanges"); biocLite("rtracklayer");'

## version of Gviz modified to show time instead of genomics units 
RUN R -e  'devtools::install_github("JoseEspinosa/Gviz")'
