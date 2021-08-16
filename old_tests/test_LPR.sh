#!/bin/bash

docker build . -t srp-analytics

docker run -v /tmp:$PWD srp-analytics test_input/7_PAH_zf_morphology_data_2020NOV11_tall.csv --LPR test_input/7_PAH_zf_LPR_data_2021JAN11_tall.csv --devel
