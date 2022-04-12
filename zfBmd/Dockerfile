FROM python:3.7
# To setup other dependencies
RUN apt-get update && apt-get install -y net-tools \
        curl \
        unixodbc \
        unixodbc-dev

COPY requirements.txt /requirements.txt
RUN pip3 install -r /requirements.txt

COPY . zfBmd
WORKDIR zfBmd

ENTRYPOINT ["python3", "main.py"]
VOLUME ["/tmp"]
