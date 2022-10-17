FROM rocker/tidyverse
# To setup other dependencies

#COPY setup.sh /setup.sh
RUN Rscript -e "install.packages('argparse')"

COPY . bmd2Samps
WORKDIR bmd2Samps

ENTRYPOINT ["Rscript","buildv2database.R"]
VOLUME ["/tmp"]
