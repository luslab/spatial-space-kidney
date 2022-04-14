FROM  rocker/tidyverse:4.0.4

# Install apt packages
RUN apt-get update \
      && apt-get install -y --no-install-recommends \
      build-essential \
      libbz2-dev \
      libz-dev
    #   apt-utils \
    #   libhdf5-dev \
    #   libboost-filesystem1.67-dev \
    #   libxtst6 \
    #   libxt-dev \
    #   libfontconfig1-dev \
    #   libfribidi-dev

# Add conda to $PATH
ENV PATH /opt/conda/bin:$PATH

# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_4.9.2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Update conda
RUN conda update -n base -c defaults conda

# Install conda packages
COPY ./environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/packages/bin:$PATH

# Setup R env
RUN head -n -1 /etc/rstudio/rserver.conf > /etc/rstudio/rserver.conf
RUN echo 'rsession-which-r=/opt/conda/envs/packages/bin/R' >> /etc/rstudio/rserver.conf

# Install R packages
RUN conda install --name packages -c conda-forge r-base
RUN conda install --name packages -c conda-forge r-devtools
RUN conda install --name packages -c conda-forge r-biocmanager
RUN conda install --name packages -c conda-forge r-tidyverse
RUN conda install --name packages -c conda-forge r-dplyr
RUN conda install --name packages -c conda-forge r-formatr
RUN conda install --name packages -c conda-forge r-remotes
RUN conda install --name packages -c conda-forge r-catools
RUN conda install --name packages -c conda-forge r-seurat
RUN conda install --name packages -c conda-forge r-patchwork

# Install Bioconductor packages
RUN R -e 'BiocManager::install("EnsDb.Mmusculus.v79")'
RUN R -e 'BiocManager::install("EnsDb.Hsapiens.v86")'
RUN R -e 'BiocManager::install("biomaRt")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("apeglm")'

# Install git builds
RUN R -e "devtools::install_github('dmcable/spacexr', build_vignettes = FALSE)"



# ENV R_HOME /opt/conda/envs/packages/bin/R

# Set default python for reticulate
# RUN   echo 'RETICULATE_PYTHON = /opt/conda/bin/python' >> /usr/local/lib/R/etc/Renviron

# Set CRAN and BioC repositories for package installations
# RUN echo 'local({library(BiocManager); r <- repositories(); r["CRAN"] <- "https://cran.rstudio.com/"; options(repos = r)}); \
#      options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))' >> /usr/local/lib/R/etc/Rprofile.site

# RUN R -e 'r = getOption("repos"); r["CRAN"] <- "https://cran.rstudio.com/"; options(repos = r)'

# Install R Packages
# RUN R -e "install.packages('devtools', repos = 'https://cran.rstudio.com/')"
    #   R -e "install.packages('future')" && \
    #   R -e "install.packages('cowplot')"  && \
#       R -e "install.packages('clustree')" && \
#       R -e "install.packages('gridExtra')" && \
#       R -e "install.packages('getopt')" && \
#       R -e "install.packages('optparse')" && \
#       R -e "install.packages('pheatmap')" && \
#       R -e "install.packages('Seurat')" && \
#       R -e "install.packages('qdapTools')" && \
#       R -e "install.packages('patchwork')" && \
#       R -e "install.packages('ggnewscale')" && \
#       R -e "remotes::install_github('mojaveazure/seurat-disk')" && \
#       # R -e "remotes::install_github('satijalab/seurat-wrappers', upgrade = 'default')" && \
#       R -e "BiocManager::install('limma')" && \
#       R -e "remotes::install_github('carmonalab/STACAS@1.1.0')" && \
#       R -e "BiocManager::install('ComplexHeatmap')" && \
#       R -e "install.packages('patchwork')" && \
#       R -e "install.packages('ggnewscale')"

# Install github R packages
# RUN R -e "devtools::install_github('dmcable/spacexr', build_vignettes = FALSE)"
# RUN   R -e "devtools::install_version('ArrayBin', version = '0.2')" && \
#       R -e "devtools::install_version('heatmap.plus', version = '1.3')" && \
#       R -e "devtools::install_github('juliendelile/Antler', dependencies = TRUE)"


