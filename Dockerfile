from rocker/tidyverse

RUN apt-get install -y net-tools\\\\\\\\\\\\\\
RUN apt-get update -qq
##        && apt-get install -y wget

COPY . srpAnalytics
WORKDIR srpAnalytics




