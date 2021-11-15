#!/bin/bash

#build docker image here!
build="/usr/bin/docker build . -t srp-analytics"
echo $build
$build

mkdir temp

##phase 3 data
#https://www.dropbox.com/sh/zg0q6wl13a3uo99/AAA0cdAK_fJwkJqpvF_HH6DWa?dl=0/

p3_morph='wget https://www.dropbox.com/sh/zg0q6wl13a3uo99/AACFZprOKkbvydjfoDI3oZo-a/Tanguay%20Phase%203%20zf%20morphology%20data%20PNNL%202021MAR23.csv -O temp/morph1.csv'
p3_lpr='wget https://www.dropbox.com/sh/zg0q6wl13a3uo99/AADd1QBStMguW9qYgzH2eatJa/Tanguay%20Phase%203%20zf%20LPR%20data%20PNNL%202021MAR23.csv -O temp/lpr1.csv'

echo $p3_morph
$p3_morph

echo $p3_lpr
$p3_lpr

##downloading files from Lisa's dropbox PFAS data
#epr_cmd="wget https://www.dropbox.com/sh/69ootcq7yyvvx2h/AABgnmHtboM4LevxK1yxPIK-a/zf%20EPA%20PFAS%20EPR_PNNL_05-28-2021.csv" -O temp/epr.csv
lpr_cmd='wget https://www.dropbox.com/sh/69ootcq7yyvvx2h/AABgzjaRPteU1EZIhnW9zv2Ka/zf%20EPA%20PFAS%20LPR_PNNL_05-28-2021.csv -O temp/lpr2.csv'
mor_cmd='wget https://www.dropbox.com/sh/69ootcq7yyvvx2h/AABxsOLgwlv7-_HTZ0xaAIlNa/zf%20EPA%20PFAS%20morphology_PNNL_05-28-2021.csv -O temp/morph2.csv'


echo $lpr_cmd
#$lpr_cmd

echo $mor_cmd
#$mor_cmd

##now run the docker image with the files
echo "Running all code"
run="/usr/bin/docker run -v "$PWD":/tmp srp-analytics --morpho=/tmp/temp/morph1.csv,/tmp/temp/morph2.csv --LPR=/tmp/temp/lpr1.csv,/tmp/temp/lpr2.csv"
echo $run
$run
