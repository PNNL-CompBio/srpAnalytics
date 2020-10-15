from rocker/tidyverse

RUN apt-get update -qq && apt-get install -y net-tools

COPY . srpAnalytics
WORKDIR srpAnalytics


CMD ["Rscript","fromEndpointsToDataFiles.R"]

VOLUME ["/tmp"]
