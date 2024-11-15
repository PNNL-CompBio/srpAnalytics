FROM python:3.10

RUN apt-get update && apt-get install -y net-tools \
        curl \
        unixodbc \
        unixodbc-dev

COPY requirements.txt /requirements.txt
RUN /usr/local/bin/python -m pip install --upgrade pip
RUN pip3 install -r /requirements.txt

COPY ./main.py zfBmd/main.py
COPY ./support_functions.py zfBMD/support_functions.py
WORKDIR zfBmd

ENTRYPOINT ["python3", "main.py"]
VOLUME ["/tmp"]
