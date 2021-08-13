# Superfund Research Program Analytics

<img src="OSU-PNNLsuperfund_Small.png"  width="400">
This repository contains the code necessary to process any new data for the Superfund Research Program Analytics Portal. Currently the portal displays two types of data:
- Zebrafish measurements describing the response to zebrafish under various levels of chemical stressors
- Environmental sample measurements that describe the relative concentration of specific chemicals in environmental samples.

Furthermore we are in the process of adding two more types of data:
- Human differential expression measurements for each chemical
- Human wristband measurements of chemical concentrations.

The data can be browsed at http://srp.pnnl.gov

This repository contains the code to handle various aspects of this portal, each described below.

## Docker image
All code in this repository requires specific package components that are contained in a Docker image. The Docker image is built automatically and stored [on Docker Hub](https://hub.docker.com/repository/docker/sgosline/srp-analytics). As such it can be pulled locally using the following command:

``` bash
docker pull sgosline/srp-analytics
```

If you would like to build the docker image locally, you can check out this repository and run build it using the following command:

``` bash
docker build . -t sgosline/srp-analytics
```


## Benchmark Dose Calculation
Calculating the benchmark dose of each chemical on the zebrafish is an active area of research. This analysis is described in an upcoming manuscript and is primarily contained in the [qc_BMD](./qc_BMD) directory. The data format required as input to this is described in the [processing pipeline schema](./schemas/processingPipelineSchema.xlsx).

Any changes to the BMD calculation will have to pass a series of tests to ensure that they work with the existing format data format


TODO: add single test

## Linking zebrafish data to environmental sample data

##
Incoming Data
We accept two types of data, each with a slightly different schema. For both types of data there is a chemical identifier (`chemical.id` or `Chemical_ID`) that must be mapped to the CAS identifiers. The schema for each type of data (the environmental sample data and the zebrafish data) can be found in the  [dataSchema spreadsheet](./dataSchemas/processingPipelineSchema.xlsx) file.

We have currently started to run this with test data to evaluate its performance.

## To run the pipeline
We have combined the data standardization to run as a single script within a docker image. This can be downloaded from Docker hub or built locally.
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

To validate an output CSV with a schema, use the following format: ```python3 validate.py <path to CSV file> <schema>```

Allowable schemas: chemdoseResponseVals, chemicalsByExtractSample, chemSummaryStats, chemXYcoords, envSampdoseResponseVals, envSampSummaryStats, or envSampXYcoords
Examples:

* ```python3 validate.py out/chemXYcoords.csv chemXYcoords```
* ```python3 validate.py out/chemdoseResponseVals.csv chemdoseResponseVals```

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
