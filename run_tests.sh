#!/bin/bash
#build database script
build="/usr/bin/docker build . -t srp-analytics"
echo $build
$build

echo "Running db test"
run="/usr/bin/docker run -v "$PWD":/tmp srp-analytics"
echo $run
$run

echo "Running validate test"
run = "/usr/bin/docker run -v "$PWD":/tmp srp-analytics --validate"
echo $run
$run

echo "Running morpho test"
run = "/usr/bin/docker run -v "$PWD":/tmp srp-analytics --test-morpho"
echo $run
$run


echo "Running LPR test"
run = '/usr/bin/docker run -v '$PWD":/tmp srp-analytics --test-lpr"
echo $run
$run
