FROM python:3.9

RUN apt-get update && apt-get install -y net-tools \
        curl \
        unixodbc \
        unixodbc-dev

COPY requirements.txt /requirements.txt
RUN /usr/local/bin/python -m pip install --upgrade pip
RUN pip3 install -r /requirements.txt

ADD *.py .
ADD srp_build_files.csv .

ENTRYPOINT ["python3", "get_build_files.py"]
VOLUME ["/tmp"]
