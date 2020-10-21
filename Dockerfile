from rocker/tidyverse

RUN apt-get update -qq && apt-get install -y net-tools
RUN Rscript -e "install.packages('argparse'); install.packages('WikipediR')"

COPY . srpAnalytics
WORKDIR srpAnalytics


CMD ["Rscript","fromEndpointsToDataFiles.R"]

VOLUME ["/tmp"]
