FROM r-base:4.2.1
# To setup other dependencies

#COPY setup.sh /setup.sh
RUN R -e "install.packages(c('argparse','tidyverse'), dependencies=TRUE,repos='http://cran.rstudio.com');"

COPY . bmd2Samps
WORKDIR bmd2Samps

ENTRYPOINT ["Rscript","buildv3database.R"]
VOLUME ["/tmp"]
