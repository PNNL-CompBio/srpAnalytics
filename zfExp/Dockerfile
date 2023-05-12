FROM bioconductor/bioconductor_docker
# To setup other dependencies

#COPY setup.sh /setup.sh
RUN Rscript -e "install.packages('readxl',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages('tidyr',dependencies=TRUE, repos='http://cran.rstudio.com')"

COPY . zfExp
COPY chemicals.csv zfExp
WORKDIR zfExp

ENTRYPOINT ["Rscript","parseGexData.R"]
VOLUME ["/tmp"]
