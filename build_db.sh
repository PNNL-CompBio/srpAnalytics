#!/bin/bash
#build database script
build="/usr/bin/docker build . -t srp-analytics"
echo $build
$build
f2=to_be_processed/7_PAH_zf_morphology_data_2020NOV11_tall.csv
f1=to_be_processed/7_PAH_zf_morphology_data_2021_Jan11.csv#
run="/usr/bin/docker run -v "$PWD":/tmp srp-analytics "$f2","$f1

echo $run
$run
