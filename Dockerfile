FROM rocker/tidyverse

COPY setup.sh setup.sh
RUN sh setup.sh

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . srpAnalytics
WORKDIR srpAnalytics

ENTRYPOINT ["python3", "dataQcBmd.py"]
VOLUME ["/tmp"]
