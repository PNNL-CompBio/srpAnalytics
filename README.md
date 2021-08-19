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

### Docker image testing
Currently all changes to the repository trigger a build of the docker image and pushing it to DockerHub. If this fails you will be notified.

If you would like to build the docker image locally, you can check out this repository and run build it using the following command:

``` bash
docker build . -t sgosline/srp-analytics
```

This is required for local testing of the code.

## Benchmark Dose Calculation
Calculating the benchmark dose of each chemical on the zebrafish is an active area of research. This analysis is described in an upcoming manuscript and is primarily contained in the [qc_BMD](./qc_BMD) directory. The data format required as input to this is described in the [processing pipeline schema](./schemas/processingPipelineSchema.xlsx).

Any changes to the BMD calculation will have to pass a series of tests to ensure that they work with the existing format data format.

### BMD Testing
Currently there are two tests for the BMD calculation, one for the morpohological changes and one for the light response. These are both automated in the continuous integration tests, but can be evaluated locally using the following commands:

``` bash
python3 dataQcBmd.py --test-morpho
python3 dataQcBmd.py --test-lpr
```

Note: these commands currently only work for chemical dose-response values, and still need to be updated to work with extract dose-response values.

## Linking zebrafish data to environmental sample data

Once we have re-calculated BMD values, we must link these data to environmental sample data. This data is stored locally on this repostiory in the [data](./data) directory, yet is read in and harmonized with the chemical data for final consumption.

### Incoming Data
The environmental sample data has a very specific format that is defined in our [schemas](./schemas) directory. To test the building the database with new data you can simply run:

``` bash
python3 dataQcBmd.py
```
This is currently being run upon pushing changes to the repository.

### Outgoing Data

To validate an output CSV with a schema, use the following format:
```
python3 validate.py <path to CSV file> <schema>
```

Allowable schemas: chemdoseResponseVals, chemicalsByExtractSample, chemSummaryStats, chemXYcoords, envSampdoseResponseVals, envSampSummaryStats, or envSampXYcoords
Examples:

* ```python3 validate.py out/chemXYcoords.csv chemXYcoords```
* ```python3 validate.py out/chemdoseResponseVals.csv chemdoseResponseVals```
