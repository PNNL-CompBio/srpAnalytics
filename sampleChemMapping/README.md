## BMD to Sample mapping

This module consumes bench mark dose measurements from the pipeline and maps them to specific chemical metadata and concentrations in sample extracts.

To run/test in this repo, navigate to the root of the repo and run:

```bash
python build_script.py --samps
```

To run/test docker image (from root of repo):

```bash
docker build . -t srp-samplechem -f sampleChemMapping/Dockerfile
docker run -v $PWD/tmp:/tmp srp-samplechem --output_dir /tmp
```

By default, this writes all output files to a `tmp` directory inside the current working directory from which the container is run. You can change this directory to point to anywhere you want output files to be written:

```bash
docker run -v /path/to/output_dir:/tmp srp-samplechem --output_dir /tmp
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
