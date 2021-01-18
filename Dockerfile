from rocker/tidyverse

RUN Rscript -e "install.packages('argparse'); install.packages('WikipediR')"
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
CMD ["python3.7","dataQcBmd.py"]

VOLUME ["/tmp"]
