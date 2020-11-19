# Superfund Research Program Analytics
<img src="OSU-PNNLsuperfund_Small.png"  width="400">
This repository contains the code necessary to process any new data for the Superfund Research Program Analytics Portal.

## Incoming Data
The incoming data will be formatted as a table with schema TBD. Currently we are just compiling the current data and not adding to it

## How to build your own
If you prefer to build your own docker image, run these commands.

```
git clone 
docker build . -t srp-analytics
docker run --volume $PWD:/tmp -ti srp-analytics
```


## How to run from docker hub
To create a data package, you simply need to add your data to the existing repository by running the following

```
docker pull sgosline/srp-analytics
```

```
docker run --volume $PWD:/tmp -ti sgosline/srp-analytics
```

The results will be the four files for the data portal.

## Data Output
The result of the pipeline. 
