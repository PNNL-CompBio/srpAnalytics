FROM python:3.9

RUN apt-get update && apt-get install -y net-tools \
        curl \
        unixodbc \
        unixodbc-dev

COPY requirements.txt /requirements.txt
RUN /usr/local/bin/python -m pip install --upgrade pip
RUN pip3 install -r /requirements.txt

COPY . zfBmd
WORKDIR zfBmd

ENTRYPOINT ["python3", "main.py"]
VOLUME ["/tmp"]