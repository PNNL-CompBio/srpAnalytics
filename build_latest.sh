#!/bin/bash

#build docker image here!
build="/usr/bin/docker build . -t srp-analytics"
echo $build
$build

mkdir temp

##downloading files from Lisa's dropbox
#epr_cmd="wget https://www.dropbox.com/sh/69ootcq7yyvvx2h/AABgnmHtboM4LevxK1yxPIK-a/zf%20EPA%20PFAS%20EPR_PNNL_05-28-2021.csv" -O temp/epr.csv
lpr_cmd='wget https://www.dropbox.com/sh/69ootcq7yyvvx2h/AABgzjaRPteU1EZIhnW9zv2Ka/zf%20EPA%20PFAS%20LPR_PNNL_05-28-2021.csv -O temp/lpr.csv'
mor_cmd='wget https://www.dropbox.com/sh/69ootcq7yyvvx2h/AABxsOLgwlv7-_HTZ0xaAIlNa/zf%20EPA%20PFAS%20morphology_PNNL_05-28-2021.csv -O temp/morph.csv'

#echo $epr_cmd
#$epr_cmd

echo $lpr_cmd
$lpr_cmd

echo $mor_cmd
$mor_cmd

##now run the docker image with the files
echo "Running all code"
run="/usr/bin/docker run -v "$PWD":/tmp srp-analytics --morpho=/tmp/temp/morph.csv --LPR=/tmp/temp/lpr.csv --get-genes"
echo $run
$run
