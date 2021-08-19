#!/bin/bash
#build database script

docker_FILE=/usr/bin/docker # Sara seems to use
if test -f "$docker_FILE"; then
    echo "$docker_FILE exists."
else 
    echo "$docker_FILE does not exist."
    docker_FILE=/usr/local/bin/docker # Doo Nam uses
fi


build=$docker_FILE" build . -t srp-analytics"
echo $build
$build

echo "Running db test"
run=$docker_FILE" run -v "$PWD":/tmp srp-analytics"
echo $run
$run

echo "Running validate test"
run=$docker_FILE" run -v "$PWD":/tmp srp-analytics --validate"
echo $run
$run

echo "Running morpho test"
run=$docker_FILE" run -v "$PWD":/tmp srp-analytics --test-morpho"
echo $run
$run


echo "Running LPR test"
run=$docker_FILE" run -v "$PWD":/tmp srp-analytics --test-lpr"
echo $run
$run
