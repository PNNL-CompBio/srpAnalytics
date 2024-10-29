## BMD to Sample mapping

This module consumes bench mark dose measurements from the pipeline and maps them to specific chemical metadata and concentrations in sample extracts. 

To run/test in this repo:
```
Rscript sampleChemMapping/mapSamplesToChems.R --sampId="+smap+' --chemId='+cid+\
            ' --epMap='+emap+' --chemClass='+cclass+\
            ' --compToxFile='+ctfile+' --sampleFiles='+fses+' --chemDesc='+desfile+\
            ' --sampMap='+smap
```

To run/test docker image (from root of repo):

```
docker build . -t srp-samplechem -f sampleChemMapping/Dockerfile
docker run -v $PWD:/tmp srp-samplechem
```

### Chemical identifiers

All chemicals must have a `cas_number` that has data download from the [EPA Comptox website](https://comptox.epa.gov/dashboard/batch-search). 

All cas numbers must have a `Chemical_ID` - these either come from the Tanguay lab or are generated manually by us. Either way they are stored in [./data/chemicalIdMapping.csv](data/chemicalIdMapping.csv). 

### Sample measurements

All sample measurements must comply with our pre-determined sample
schema. There must be a `Sample_ID` mapping to `SampleNumber` in the
sample mapping file.

### Benchmark dose values

These are processed from the stored BMD files and the recalculated
ones passed into the argument. 
