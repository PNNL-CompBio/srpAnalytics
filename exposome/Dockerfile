FROM r-base:4.4.1
# To setup other dependencies


RUN apt-get update --allow-insecure-repositories
RUN apt-get install -y --allow-unauthenticated --fix-missing  python3-pip python3-setuptools python3-dev python3-venv libcurl4-openssl-dev libglpk-dev libxml2-dev libpq-dev


RUN python3 -m venv /opt/venv

ENV PATH="/opt/venv/bin:$PATH"
ENV PYTHONPATH "${PYTHONPATH}:/tmp"
ENV PYTHONPATH "${PYTHONPATH}:/app"
WORKDIR /app

RUN mkdir -p exposome
COPY exposome/* exposome/
COPY build_script.py ./
COPY srpAnalytics.yaml ./
COPY requirements.txt .

RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt
RUN Rscript -e "install.packages(c('jsonlite','dplyr','httr','tidyr'), dependencies=TRUE,repos='http://cran.rstudio.com')"


ENTRYPOINT ["python3","/app/build_script.py","--expo"]
