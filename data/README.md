## Data directory
This directory contains existing data and references files to be added to the SRP analytics database. These data are not private in any way, but mainly serve as reference files for the broader repository.

### Data file manifest
The data file manifest contains a link to all the other data files,
including those referenced in this document. This file can be found
[here](./srp_build_files.csv). 


### Data Schema and validaton
Incoming data requires a pre-defined schema so we can validate things
BEFORE they are going to be processed by our pipeline. This schema
also enables an easy way to check if a file CAN be processed before
adding it to the list.

The data schema file can be found [here](./srp_templates.yaml). 

To validate incoming data your file must be one of the following
types:

1. morphology
2. behavior
3. zebrafishMapping
4. fses

You must download install the linkML tool using `pip` by following the
instructions [on the LinkML
page](https://linkml.io/linkml/intro/install.html).

To validate a single file you need the type of file and the path to
the file.

```
linkml-validate fses/FSES_indoor_outdoor_study.csv -C fses -s srp_templates.yaml
```

This will validate the fses file.

### Chemical reference information

Many of the build files are actual reference information. 

#### Chemical id file
#### Chemical classification file
#### compTox file

### Environmental sample information

#### Sample id file

### Existing Zebrafish data


####
