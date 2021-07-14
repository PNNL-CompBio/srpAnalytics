# Superfund Research Program Analytics
<img src="OSU-PNNLsuperfund_Small.png"  width="400">
This repository contains the code necessary to process any new data for the Superfund Research Program Analytics Portal. It is designed to consume two types of data: (1) zebrafish measurements describing the response to zebrafish under various levels of chemical stressors, and (2) environmental sample measurements that describe the relative concentration of specific chemicals in environmental samples. This is described in more detail in our `Scientific Data` manuscript (to be submitted) and depicted below.
<img src='overview.jpg' width="200">

## Incoming Data
We accept two types of data, each with a slightly different schema. For both types of data there is a chemical identifier (`chemical.id` or `Chemical_ID`) that must be mapped to the CAS identifiers. The schema for each type of data (the environmental sample data and the zebrafish data) can be found in the  [dataSchema spreadsheet](./dataSchemas/processingPipelineSchema.xlsx) file.

We have currently started to run this with test data to evaluate its performance.

## To run the pipeline
We have combined the data standardization
If you prefer to build your own docker image, run these commands.

```
git clone https://github.com/sgosline/srpAnalytics.git
cd srp-analytics
docker build . -t srp-analytics
```
(this building took 8 minutes in mackbook since it installs all depencies)

Then run like this
```
docker run -v $PWD:/tmp srp-analytics [your file here, for example to_be_processedt/7_PAH_zf_morphology_data_2020NOV11_tall.csv]
```

To build the whole database:

``` 1c-enterprise
sh build_db.sh
```

## How to run from docker hub
To create a data package, you simply need to add your data to the existing repository by running the following

```
docker pull sgosline/srp-analytics
```
(this pulling took 5 minutes in mackbook

Then run like this
```
docker run -v $PWD:/tmp sgosline/srp-analytics [your file here, for example test_input/7_PAH_zf_morphology_data_2020NOV11_tall.csv]
```

The results will be the four files for the data portal. Add the `--devel` flag if you are just testing the code.
