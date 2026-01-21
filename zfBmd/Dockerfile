FROM python:3.10

RUN apt-get update && apt-get install -y net-tools \
        curl \
        unixodbc \
        unixodbc-dev

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
WORKDIR /app

RUN mkdir ./zfBmd
COPY zfBmd ./zfBmd
RUN mkdir data
COPY data/srp_build_files.csv data/srp_build_files.csv
COPY srpAnalytics.yaml .
COPY zfBmd/requirements.txt .
COPY build_script.py .
RUN mkdir -p ./tmp/.figshare_cache

RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt

ENTRYPOINT ["python3", "-u", "/app/build_script.py", "--bmd"]
