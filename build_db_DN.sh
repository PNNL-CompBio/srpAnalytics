#!/bin/bash
#build database script
build="docker build . -t srp-analytics"
echo $build
$build

#f1=to_be_processed/7_PAH_zf_morphology_data_2020NOV11_tall.csv
f1=to_be_processed/7_PAH_zf_morphology_data_2020NOV11_tall_3756.csv

run="docker run -v "$PWD":/tmp srp-analytics "$f1" --devel"
echo $run
$run
