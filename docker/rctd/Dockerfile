FROM jupyter/r-notebook

# # Install conda packages
# COPY ./environment.yml /
# RUN conda env create -f /environment.yml && conda clean -a
# ENV PATH /opt/conda/envs/reporting/bin:$PATH

# Install R Packages
RUN R -e "devtools::install_github('dmcable/spacexr', build_vignettes = FALSE)"