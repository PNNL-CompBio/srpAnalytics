FROM rocker/tidyverse
# To setup other dependencies


RUN apt-get update --allow-insecure-repositories
RUN apt-get install -y --allow-unauthenticated --fix-missing gcc-12 build-essential libpq-dev python3-venv


RUN python3 -m venv /opt/venv


ENV PATH="/opt/venv/bin:$PATH"
ENV PYTHONPATH "${PYTHONPATH}:/tmp"
ENV PYTHONPATH "${PYTHONPATH}:/app"
WORKDIR /app

RUN mkdir -p sapleChemMapping
COPY sampleChemMapping/* sampleChemMapping
COPY build_script.py .
COPY requirements.txt .

RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt

RUN R -e "install.packages(c('argparse','tidyr','dplyr','rio','xml2'), dependencies=TRUE,repos='http://cran.rstudio.com');"


ENTRYPOINT ["python","/app/build_script.py","--samps"]

