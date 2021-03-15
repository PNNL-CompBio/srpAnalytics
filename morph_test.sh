#!/bin/bash
docker build . -t srp-analytics
docker run -v $PWD:/tmp srp-analytics test_input/7_PAH_zf_morphology_data_2020NOV11_tall.csv --devel
