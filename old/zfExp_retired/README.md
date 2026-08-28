## ZF Expression parsing

This module parses the ZF expression data into files amenable to the
data portal.

This module also runs the `build_script.py` script to build the zebrafish
expression data.

```bash
python build_script.py --geneEx
```

However the build script is also copied into the Dockerfile and run as below:

```bash
docker build . -t srp-zfexp -f zfExp/Dockerfile
docker run -v $PWD/tmp:/tmp srp-zfexp
```

By default, this writes all output files to a `tmp` directory inside the current working directory from which the container is run. You can change this directory to point to anywhere you want output files to be written:

```bash
docker run -v /path/to/output_dir:/tmp srp-zfexp --output_dir /tmp
```
