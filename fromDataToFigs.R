##now that we have behavioral data and Serberus platform this script is designed to
#process everything into a digestable data

require(dplyr)
require(readxl)

summ = read.csv('chemSummaryStats.csv')
#read.csv(curves,file='/tmp/chemXYcoords.csv',row.names = FALSE)
#read.csv(doseReps,file='/tmp/chemdoseResponseVals.csv',row.names = FALSE)
#write.csv(sampChem,file='/tmp/chemicalsByExtractSample.csv',row.names=FALSE)
ext <-read.csv('chemicalsByExtractSample.csv')

chemsOnly<-summ%>%left_join(ext,by='Chemical_ID')%>%distinct()%>%
    #subset(Model!='NULL')%>%
   # subset(measurement_value>0.0)%>%
    group_by(chemical_class)%>%
    summarize(`Number of Chemicals`=n_distinct(Chemical_ID),
              `Number of endpoints`=n_distinct(End_Point[which(Model!='NULL')]),
              `Number of extracts`=n_distinct(extract_exp[which(measurement_value>0.0)]))

##do we have curves for the extract data? if so we need to merge that as well
eps<-summ%>%select(End_Point,End.Point.Name,Description,endPointLink)%>%distinct()

write.csv(chemsOnly,'chemCounts.csv',row.names = F)
write.csv(eps,'epDesc.csv',row.names=F)