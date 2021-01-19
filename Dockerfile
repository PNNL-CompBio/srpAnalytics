<<<<<<< HEAD
from rocker/tidyverse

RUN apt-get update -qq && apt-get install -y net-tools
RUN Rscript -e "install.packages('argparse'); install.packages('WikipediR')"

COPY . srpAnalytics
WORKDIR srpAnalytics


CMD ["Rscript","fromEndpointsToDataFiles.R"]

VOLUME ["/tmp"]
=======
from rocker/tidyverse

RUN Rscript -e "install.packages('argparse')"
RUN apt-get update -qq && apt-get install -y net-tools
RUN apt-get install -y python3.7
RUN apt-get install -y python3-pip
RUN pip3 install pandas
RUN pip3 install matplotlib
RUN pip3 install seaborn
RUN pip3 install statsmodels
RUN pip3 install astropy
RUN pip3 install scipy==1.4.1

COPY . srpAnalytics
WORKDIR srpAnalytics

#CMD ["Rscript","fromEndpointsToDataFiles.R"]
#CMD ["python3","dataQcBmd.py"]

VOLUME ["/tmp"]
>>>>>>> 93d5103ae653005638b9126e9865ec11c653de1a
