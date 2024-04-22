## BMD to Sample mapping

This module consumes bench mark dose measurements from the pipeline and maps them to specific chemical metadata and concentrations in sample extracts. 

### Chemical identifiers

All chemicals must have a `cas_number` that has data download from the [EPA Comptox website](https://comptox.epa.gov/dashboard/batch-search). 

All cas numbers must have a `Chemical_ID` - these either come from the Tanguay lab or are generated manually by us. Either way they are stored in [./data/chemicalIdMapping.csv](data/chemicalIdMapping.csv). 

### Sample measurements

All sample measurements must comply with our pre-determined sample schema.

### Benchmark dose values