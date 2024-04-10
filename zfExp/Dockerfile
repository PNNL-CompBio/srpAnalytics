FROM bioconductor/bioconductor_docker
# To setup other dependencies

#COPY setup.sh /setup.sh
RUN Rscript -e "install.packages('readxl',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages('tidyr',dependencies=TRUE, repos='http://cran.rstudio.com')"

RUN mkdir -p zfExp
COPY parseGexData.R zfExp/
#COPY chemicals.csv zfExp
WORKDIR zfExp

ENTRYPOINT ["Rscript","parseGexData.R"]
VOLUME ["/tmp"]
