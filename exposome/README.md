## Exposome Data Processing

This module contains the scripts and docker image to build the
exposome data.

To test locally, you can run this script:

```bash
python build_script.py --expo

```

To test with all the R/Python requirements already installed, use the
docker image:

```bash
docker build . -t srp-exposome -f exposome/Dockerfile
docker run -v $PWD/tmp:/tmp srp-exposome
```

By default, this writes all output files to a `tmp` directory inside the current working directory from which the container is run. You can change this directory to point to anywhere you want output files to be written:

```bash
docker run -v /path/to/output_dir:/tmp srp-exposome --output_dir /tmp
```
