# Benchmark Dose-Response Curve (BMDRC) Fitting

To build the Docker container, navigate to the top level directory of the repository (i.e. from this directory, run `cd ..`) and run:

```bash
docker build . -t srp-zfbmd -f zfBmd/Dockerfile
```

To run the container on the provided test data, use:

```bash
# Morphology example
docker run -v $PWD/tmp:/tmp srp-zfbmd --morpho /app/zfBmd/test_files/test_morphology.csv --output /tmp

# LPR example
docker run -v $PWD/tmp:/tmp srp-zfbmd --lpr /app/zfBmd/test_files/test_behavioral.csv --output /tmp

# Both
docker run -v $PWD/tmp:/tmp srp-zfbmd --morpho /app/zfBmd/test_files/test_morphology.csv --lpr /app/zfBmd/test_files/test_behavioral.csv --output /tmp
```

By default, this writes all output files to a `tmp` directory inside the current working directory from which the container is run. You can change this directory to point to anywhere you want output files to be written:

```bash
docker run -v /path/to/output_dir:/tmp srp-zfbmd ...FLAGS... --output /tmp
```
