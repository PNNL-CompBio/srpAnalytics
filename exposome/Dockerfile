FROM rocker/tidyverse
# To setup other dependencies

#COPY setup.sh /setup.sh
RUN Rscript -e "install.packages('argparse')"

COPY . exposome
WORKDIR exposome

ENTRYPOINT ["Rscript","exposome_summary_stats.R"]
VOLUME ["/tmp"]
