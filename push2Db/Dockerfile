FROM python:3.7
# To setup other dependencies
RUN apt-get update -qq && apt-get install -y net-tools \
        curl \
        unixodbc \
        unixodbc-dev

COPY requirements.txt /requirements.txt
RUN pip3 install -r /requirements.txt

COPY . push2Db
WORKDIR push2Db

ENTRYPOINT ["python3", "main.py"]
VOLUME ["/tmp"]
