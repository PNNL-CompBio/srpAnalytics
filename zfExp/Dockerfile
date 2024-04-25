FROM r-base:4.2.1
# To setup other dependencies


RUN apt-get update --allow-insecure-repositories
RUN apt-get install -y --allow-unauthenticated --fix-missing gcc-12 build-essential libpq-dev python3-venv

RUN python3 -m venv /opt/venv

ENV PATH="/opt/venv/bin:$PATH"
ENV PYTHONPATH "${PYTHONPATH}:/tmp"
ENV PYTHONPATH "${PYTHONPATH}:/app"

WORKDIR /app
RUN mkdir -p zfExp
COPY zfExp/* zfExp/
COPY build_script.py ./
COPY requirements.txt ./

RUN Rscript -e "install.packages('readxl',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages('tidyr',dependencies=TRUE, repos='http://cran.rstudio.com')"
RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt


ENTRYPOINT ["python","/app/build_script.py","--geneEx"]

