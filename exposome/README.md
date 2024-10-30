## Exposome data processing
This module contains the scripts and docker image to build the
exposome data.


To test locally you can run this script:
```
python build_script.py --expo

```

To test with all the R/Python requirements already installed, use the
docker image:

```
docker build . -t srp-exposome -f exposome/Dockerfile
docker run -v $PWD:/tmp srp-exposome
```

    
