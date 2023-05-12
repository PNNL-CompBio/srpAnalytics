FROM bioconductor/bioconductor_docker
# To setup other dependencies

#COPY setup.sh /setup.sh
RUN R -e "install.packages('argparse', dependencies=TRUE,repos='http://cran.rstudio.com');"

COPY . bmd2Samps
WORKDIR bmd2Samps

ENTRYPOINT ["Rscript","buildv2database.R"]
VOLUME ["/tmp"]
