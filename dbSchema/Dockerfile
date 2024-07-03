FROM python:3.9
#-alpine3.18
# To setup other dependencies
RUN apt-get update -qq && apt-get install -y net-tools \
        curl \
        unixodbc \
        unixodbc-dev

COPY requirements.txt /requirements.txt
RUN pip3 install -r /requirements.txt

RUN mkdir dbSchema
COPY dbSchema/* dbSchema/
#COPY build_script.py .
WORKDIR dbSchema

ENTRYPOINT ["python3", "main.py"]
VOLUME ["/tmp"]
