# Superfund Research Program Analytics

<img src="OSU-PNNLsuperfund_Small.png"  width="400">
This repository contains the code necessary to process any new data for the Superfund Research Program Analytics Portal.

## Incoming Data

The incoming data will be formatted as a table with schema TBD. Currently we are just compiling the current data and not adding to it

## How to build your own docker image

If you prefer to build your own docker image, run these commands (turn off PNNL VPN).

```bash
git clone https://github.com/sgosline/srpAnalytics.git
cd srp-analytics
docker build . -t srp-analytics
```

(this building took 8 minutes in macbook pro since it installs all dependencies)

Then, run like this

```bash
docker run -v $PWD:/tmp srp-analytics [your file here, for example to_be_processed/7_PAH_zf_morphology_data_2020NOV11_tall.csv] [--update-db]
```

If ```--update-db``` is not included, the script will still test the connection to the database. If it is included the script, by default, overwrites the develop database. You can either append or replace to either the develop or production database (requires changes to the ```dataQcBmd.py``` file).

For faster running add --devel,

```bash
docker run -v $PWD:/tmp srp-analytics [your file here, for example to_be_processed/7_PAH_zf_morphology_data_2020NOV11_tall.csv] --devel
```

To build the whole database:

``` 1c-enterprise
sh build_db.sh
```

## How to run from docker hub

To create a data package, you simply need to add your data to the existing repository by running the following

```bash
docker pull sgosline/srp-analytics
```

(this pulling took 5 minutes in mackbook

Then run like this

```bash
docker run -v $PWD:/tmp sgosline/srp-analytics [your file here, for example test_input/7_PAH_zf_morphology_data_2020NOV11_tall.csv]
```

The results will be the four files for the data portal. Add the `--devel` flag if you are just testing the code.

## Data Input

Currently the code is designed to take a specific form of input to be processed by our benchmark dose analysis pipeline. The columns are defined in the table below:

| Column name| Description|
| ---| ---|

## Data Output

The result of the pipeline will be six files, zipped up into a single resource.
