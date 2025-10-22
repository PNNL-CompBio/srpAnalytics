# FROM r-base:4.4.1
FROM python:3.12
# To setup other dependencies

RUN apt-get update --allow-insecure-repositories
RUN apt-get install -y --allow-unauthenticated --fix-missing  python3-pip python3-setuptools python3-dev python3-venv libcurl4-openssl-dev libglpk-dev libxml2-dev libpq-dev
# r-base
#gcc-12 build-essential libpq-dev python3-venv libxml2-dev libcurl4-openssl-dev

# RUN R -e "install.packages(c('argparse','tidyr','dplyr','rio','xml2'), dependencies=TRUE,repos='http://cran.rstudio.com');"
RUN python3 -m venv /opt/venv


ENV PATH="/opt/venv/bin:$PATH"
ENV PYTHONPATH="${PYTHONPATH}:/tmp"
ENV PYTHONPATH="${PYTHONPATH}:/app"
WORKDIR /app

RUN mkdir -p sampleChemMapping
COPY sampleChemMapping/* sampleChemMapping
COPY build_script.py .
COPY srpAnalytics.yaml .
COPY requirements.txt .

RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt

ENTRYPOINT ["python","/app/build_script.py","--samps"]

