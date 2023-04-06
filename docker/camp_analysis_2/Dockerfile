FROM  rocker/tidyverse:4.0.4

# Install apt packages
RUN apt-get update \
      && apt-get install -y --no-install-recommends \
      build-essential \
      libbz2-dev \
      libz-dev

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

# Install mamba
RUN conda install --name base -c conda-forge mamba

# Install conda packages
RUN mamba create --name packages python=3.8
ENV PATH /opt/conda/envs/packages/bin:$PATH

# Setup R env
RUN head -n -1 /etc/rstudio/rserver.conf > /etc/rstudio/rserver.conf
RUN echo 'rsession-which-r=/opt/conda/envs/packages/bin/R' >> /etc/rstudio/rserver.conf

# Install R base packages
RUN mamba install --name packages -c conda-forge r-base
RUN mamba install --name packages -c conda-forge r-devtools
RUN mamba install --name packages -c conda-forge r-biocmanager
RUN mamba install --name packages -c conda-forge r-tidyverse
RUN mamba install --name packages -c conda-forge r-dplyr
RUN mamba install --name packages -c conda-forge r-formatr
RUN mamba install --name packages -c conda-forge r-remotes
RUN mamba install --name packages -c conda-forge r-catools
RUN mamba install --name packages -c conda-forge r-doparallel

# Install tool specific packages
RUN mamba install --name packages -c conda-forge r-seurat
RUN mamba install --name packages -c conda-forge r-patchwork

# Install bioconductor
# RUN mamba install --name packages -c bioconda bioconductor-ensdb.hsapiens.v86
RUN mamba install --name packages -c bioconda bioconductor-ensdb.mmusculus.v79
RUN mamba install --name packages -c bioconda bioconductor-biomart
RUN mamba install --name packages -c bioconda bioconductor-deseq2

# bioconductor-apeglm 1.20.0

# # Install Bioconductor packages
# RUN R -e 'BiocManager::install("DESeq2")'
# RUN R -e 'BiocManager::install("apeglm")'

# Install git builds
# RUN R -e "devtools::install_github('angelolab/FlowSOM', upgrade = FALSE, upgrade_dependencies = FALSE)"
RUN R -e "devtools::install_github('dmcable/spacexr', build_vignettes = FALSE)"
