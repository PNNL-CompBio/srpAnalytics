FROM python:3.10

RUN apt-get update && apt-get install -y net-tools \
        curl \
        unixodbc \
        unixodbc-dev

COPY zfBmd/requirements.txt /requirements.txt
RUN /usr/local/bin/python -m pip install --upgrade pip
RUN pip3 install -r /requirements.txt

RUN mkdir zfBmd
COPY zfBmd/main.py zfBmd/main.py
COPY zfBmd/support_functions.py zfBmd/support_functions.py
COPY zfBmd/test_files/* zfBmd/test_files/
WORKDIR zfBmd

ENTRYPOINT ["python3", "main.py"]
VOLUME ["/tmp"]
