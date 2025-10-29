FROM r-base:4.4.1
FROM python:3.12
# To setup other dependencies


RUN apt-get update --allow-insecure-repositories
RUN apt-get install -y --allow-unauthenticated --fix-missing  python3-pip python3-setuptools python3-dev python3-venv libcurl4-openssl-dev libglpk-dev libxml2-dev libpq-dev r-base
RUN apt-get install -y libmariadb-dev-compat libmariadb-dev

RUN python3 -m venv /opt/venv

ENV PATH="/opt/venv/bin:$PATH"
# ENV PYTHONPATH="${PYTHONPATH}:/tmp"
# ENV PYTHONPATH="${PYTHONPATH}:/app"

WORKDIR /app
RUN mkdir -p zfExp
COPY zfExp/* zfExp/
COPY build_script.py ./
COPY srpAnalytics.yaml ./
COPY requirements.txt ./

RUN mkdir -p sampleChemMapping
COPY sampleChemMapping/* sampleChemMapping/
RUN mkdir src
COPY src/* src

RUN Rscript -e "install.packages('readxl',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages('tidyr',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages('rio',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages('rjson',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages('WriteXLS',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/enrichR/enrichR_3.2.tar.gz', dependencies=TRUE, source=TRUE)"
RUN Rscript -e "install.packages('argparse',dependencies=TRUE,repos='http://cran.rstudio.com')"

RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt


ENTRYPOINT ["python", "/app/build_script.py", "--geneEx"]

