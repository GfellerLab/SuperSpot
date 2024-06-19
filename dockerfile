# Use the latest R version as base image
FROM rocker/tidyverse

# Set the locale environment variables
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
RUN apt-get update && apt-get install -y locales && \
    locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8

RUN apt-get update && apt-get install -y patch
RUN apt-get update && apt-get install -y gcc g++

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libhdf5-dev \
    libhdf5-serial-dev \
    hdf5-tools


# Install system dependencies for certain R packages, including GSL
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libglpk-dev \
    libgit2-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libgsl-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages from CRAN
RUN R -e "install.packages('Rhdf5lib', dependencies=TRUE)"
RUN R -e "install.packages('rhdf5filters', dependencies=TRUE)"
RUN R -e "install.packages('hdf5r', dependencies=TRUE)"


RUN R -e "install.packages(c('igraph', 'devtools','RANN', 'WeightedCluster', 'corpcor', 'weights', 'Hmisc', 'Matrix', 'matrixStats', 'plyr', 'irlba', 'grDevices', 'patchwork', 'gtools', 'ggplot2', 'umap', 'entropy', 'Rtsne', 'dbscan', 'cowplot', 'scales', 'plotfunctions', 'proxy', 'methods', 'rlang', 'Seurat', 'data.table', 'foreach', 'doParallel', 'concaveman', 'parallelDist', 'pbapply', 'pliman', 'tidyr', 'tidyverse','bench','Rfast2','hdf5r','arrow'), dependencies=TRUE)"

# Reinstall SeuratObject to match the current R version
RUN R -e "install.packages('SeuratObject', dependencies=TRUE)"

# Set repositories and install additional packages
RUN R -e "setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'));install.packages(c('BPCells', 'presto', 'glmGamPoi','peakRAM','Rfast2'))"

# Install Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('scater', 'bluster'))"

# Install GitHub packages
RUN R -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools'); devtools::install_github('GfellerLab/SuperCell')"
RUN R -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools'); devtools::install_github('GfellerLab/SuperSpot')"
RUN R -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools'); devtools::install_github('immunogenomics/presto')"

# Set default command to launch R
#CMD ["R"]
