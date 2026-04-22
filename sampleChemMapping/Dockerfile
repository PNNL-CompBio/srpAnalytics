FROM python:3.12
# To setup other dependencies

RUN apt-get update --allow-insecure-repositories
RUN apt-get install -y --allow-unauthenticated --fix-missing  python3-pip python3-setuptools python3-dev python3-venv libcurl4-openssl-dev libglpk-dev libxml2-dev libpq-dev
# r-base gcc-12 build-essential libpq-dev python3-venv libxml2-dev libcurl4-openssl-dev

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
WORKDIR /app

RUN mkdir -p sampleChemMapping
COPY sampleChemMapping/* sampleChemMapping
RUN mkdir -p src
COPY src/* src
RUN mkdir data
COPY data/srp_build_files.csv data/srp_build_files.csv
COPY build_script.py .
COPY srpAnalytics.yaml .
COPY requirements.txt .
RUN mkdir -p tmp/.figshare_cache

RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt

ENTRYPOINT ["python", "/app/build_script.py", "--samps"]

