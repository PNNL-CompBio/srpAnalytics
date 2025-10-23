# Build command: docker build -t zfbmd .
#
# Morphology example using a local folder called tmp and data in the container: 
#       docker run -v ./tmp:/zfBmd/tmp zfbmd --morpho test_files/test_morphology.csv --output /zfBmd/tmp/
#
# Both morphology and LPR example using a local folder called tmp and data in the container:
#       docker run -v ./tmp:/zfBmd/tmp zfbmd --morpho test_files/test_morphology.csv --lpr test_files/test_behavioral.csv --output /zfBmd/tmp/
#
# Both example using local data. It must be moved to the tmp folder to work
#       docker run -v ./tmp:/zfBmd/tmp zfbmd --morpho /zfBmd/tmp/test_morphology.csv --lpr /zfBmd/tmp/test_behavioral.csv --output /zfBmd/tmp/
FROM python:3.10

RUN apt-get update && apt-get install -y net-tools \
        curl \
        unixodbc \
        unixodbc-dev

COPY requirements.txt /requirements.txt
RUN /usr/local/bin/python -m pip install --upgrade pip
RUN pip3 install -r /requirements.txt

COPY ./main.py zfBmd/main.py
COPY ./support_functions.py zfBmd/support_functions.py
COPY ./test_files/* zfBmd/test_files/
WORKDIR zfBmd
RUN mkdir ./tmp

ENTRYPOINT ["python3", "-u", "main.py"]
