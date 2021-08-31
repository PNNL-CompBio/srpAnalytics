#!/bin/bash
#build database script

docker_FILE=/usr/bin/docker # Sara seems to use
if test -f "$docker_FILE"; then
    echo "$docker_FILE exists."
else 
    echo "$docker_FILE does not exist."
    docker_FILE=/usr/local/bin/docker # Doo Nam uses
fi

echo "
Building the docker image locally"
build=$docker_FILE" build . -t srp-analytics"
echo $build
$build

echo "
Running db test"
run=$docker_FILE" run -v "$PWD":/tmp srp-analytics"
echo $run
$run

echo "
Running validate test"
run=$docker_FILE" run -v "$PWD":/tmp srp-analytics --validate"
echo $run
$run


echo "
Running morpho test with 7 PAH data (single chemical)"
run=$docker_FILE" run -v "$PWD":/tmp srp-analytics --test-morpho"
echo $run
$run


echo "
Running morpho test with extract data (single chemical)"
run=$docker_FILE" run -v "$PWD":/tmp srp-analytics --test-extract"
echo $run
$run



START_TIME=$(date +%s)

echo "
Running LPR test (warning: before actual BMD calculation for this lpr test, preprocessing (e.g. reformat) takes long time like 5 minutes)"
run=$docker_FILE" run -v "$PWD":/tmp srp-analytics --test-lpr"
echo $run
$run

END_TIME=$(date +%s)

echo "
LPR test took $(($END_TIME - $START_TIME)) seconds"
# (when echoing, scp/rsync from my mac didn't work)
