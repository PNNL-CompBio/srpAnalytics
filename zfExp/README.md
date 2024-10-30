## ZF Expression parsing

This module parses the ZF expression data into files amenable to the
data portal.


This module also runs the `build_all.py` script to build the zf
expression data. 

```
python build_script.py --zfExp
```

However the build script is also copied into the Dockerfile and run as below:

```
docker build . -t srp-zfexp -f zfExp/Dockerfile
docker run -v $PWD:/tmp srp-zfexp

```
