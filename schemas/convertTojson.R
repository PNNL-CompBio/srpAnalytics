##function that converst excel to xlsx

library(readxl)
library(jsonlite)

env <- readxl::read_xlsx('processingPipelineSchema.xlsx',sheet='environmentalSample')
zeb <-readxl::read_xlsx('processingPipelineSchema.xlsx',sheet='zebrafish raw data')

write(jsonlite::toJSON(apply(env,1,as.list)),'envSampleIntake.json')
write(jsonlite::toJSON(apply(zeb,1,as.list)),'zebrafishDataIntake.json')