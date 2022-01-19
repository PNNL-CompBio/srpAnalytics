#!/bin/bash
#build database script

docker_FILE=/usr/bin/docker # Sara seems to use
if test -f "$docker_FILE"; then
    echo "$docker_FILE exists."
else
    echo "$docker_FILE does not exist."
    docker_FILE=/usr/local/bin/docker # Doo Nam uses
fi

$docker_FILE pull sgosline/srp-exposome
$docker_FILE pull sgosline/srp-zfbmd
$docker_FILE pull sgosline/srp-zf2samps

echo "Running gene and db test"
run="/usr/bin/docker run -v "$PWD":/tmp sgosline/srp-dbschema"
echo $run
#$run

echo "
Running validate test"
run=$docker_FILE" run -v "$PWD":/tmp sgosline/srp-dbschema --validate"
echo $run
#$run

echo "
Running morpho test with 7 PAH data (single chemical)"
run=$docker_FILE" run -v "$PWD":/tmp sgosline/srp-zfbmd --test-morpho"
echo $run
$run


echo "
Running morpho test with extract data (single chemical)"
run=$docker_FILE" run -v "$PWD":/tmp sgosline/srp-zfbmd --test-extract"
echo $run
$run

START_TIME=$(date +%s)

echo "
Running LPR test (warning: before actual BMD calculation for this lpr test, preprocessing (e.g. reformat) takes long time like 5 minutes)"
run=$docker_FILE" run -v "$PWD":/tmp sgosline/srp-zfbmd --test-lpr"
echo $run
$run

END_TIME=$(date +%s)

echo "
LPR test took $(($END_TIME - $START_TIME)) seconds"
# (when echoing, scp/rsync from my mac didn't work)

echo "
Combining all files into single db"
run = $docker_FILE" run -v "$PWD":/tmp sgosline/srp-bmd2samps"
echo $run
$run
