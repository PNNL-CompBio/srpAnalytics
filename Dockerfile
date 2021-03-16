from rocker/tidyverse

RUN apt-get update -qq && apt-get install -y net-tools \
        python3.7 \
        python3-pip

RUN Rscript -e "install.packages('argparse')"

RUN pip3 install pandas
RUN pip3 install matplotlib
RUN pip3 install seaborn
RUN pip3 install statsmodels
RUN pip3 install astropy
RUN pip3 install scipy==1.4.1

COPY . srpAnalytics
WORKDIR srpAnalytics

ENTRYPOINT ["python3", "dataQcBmd.py"]
VOLUME ["/tmp"]
