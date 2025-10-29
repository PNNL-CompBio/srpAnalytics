FROM python:3.10

RUN apt-get update && apt-get install -y net-tools \
        curl \
        unixodbc \
        unixodbc-dev

COPY zfBmd/requirements.txt /requirements.txt
RUN /usr/local/bin/python -m pip install --upgrade pip
RUN pip3 install -r /requirements.txt

RUN mkdir /zfBmd
COPY zfBmd /zfBmd
WORKDIR /zfBmd
RUN mkdir ./tmp

ENTRYPOINT ["python3", "-u", "main.py"]
