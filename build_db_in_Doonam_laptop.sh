#!/bin/bash
#build database script
build="/usr/local/bin/docker build . -t srp-analytics"
echo $build
$build
f2=to_be_processed/7_PAH_zf_morphology_data_2020NOV11_tall.csv

run="/usr/local/bin/docker run -v "$PWD":/tmp srp-analytics "$f2""

echo $run
$run
