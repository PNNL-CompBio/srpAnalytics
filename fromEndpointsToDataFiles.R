##now that we have behavioral data and Serberus platform this script is designed to
#process everything into a digestable data

require(dplyr)
require(readxl)
require(argparse)
require(xml2)
#require(WikipediR)
##The data release will be comprised of 5 files
#' 1- list of chemicals nad their information
#' 2- list of sites, their metadata, and the chemical composition (curated sample data)
#' 3- ZF summary statistics
#' 4- ZF points to plot
#' 5- ZF curves to plot


#These pathways refer to absolute pathways in the docker image
data.dir='/srpAnalytics/data/'
#data.dir='data/'
out.dir='/tmp/'
#out.dir='./'

#step 1 - get cas ids to pre-defined class
chemMapping<-read.csv(paste0(data.dir,'Chemicals.csv'))%>%
                      dplyr::select(cas_number,chemical_class,zaap_cid)%>%
                      distinct()%>%
                      rename(zf.cid='zaap_cid')

#this file was downloaded from osu
newMapping<-readxl::read_xlsx(paste0(data.dir,'chemicalid_table_2020AUG07.xlsx'))%>%
    select(cas_number='casrn',zf.cid='chemical_id')%>%
    mutate(chemical_class=NA)

#here is more mapping information from OSU
fullMapping<-readxl::read_xlsx(paste0(data.dir,'OSU_BU_Overlap_Inventory_2020MAR30.xlsx'))%>%
    select(casrn,zf.cid)%>%
    distinct()%>%
    tidyr::separate_rows(zf.cid,sep=';\\s+')%>%
                                        #rename(cas_number='casrn')%>%
    left_join(chemMapping)%>%
    select(-cas_number)%>%
    rename(cas_number='casrn')%>%
    rbind(newMapping)

#here is our pre-defined dictionary
endpointDetails<-readxl::read_xlsx(paste0(data.dir,'SuperEndpoint Mapping 2020Sept2.xlsx'),
                                   sheet='Dictionary')%>%
    rename(End_Point='Abbreviation',`End Point Name`='Simple name (<20char)')%>%
  rename(endPointLink='Ontology Link')


#'getNewChemicalClass
#'This file reads in the files and processes the chemical class names
#'to be friendly for the website
getNewChemicalClass<-function(){
  pahs=readxl::read_xlsx(paste0(data.dir,'/PAH_and_1530_SRP_Summary.xlsx'),sheet='Master summary-PAHs')%>%
      dplyr::select(cas_number='casrn')%>%
      mutate(Classification='PAH')
  extras<-data.frame(cas_number=c("3074-03-01","7496-02-08","6373-11-01"),Classification='PAH')
  non.pahs=readxl::read_xlsx(paste0(data.dir,'/PAH_and_1530_SRP_Summary.xlsx'),sheet='Master Summary-Non-PAHs')%>%
    dplyr::select(cas_number='casrn',Classification='classification')

  full.class<-rbind(pahs,extras,non.pahs)
  full.class$newClass=sapply(full.class$Classification,function(x)
    if(x%in%c('industrial',"industrial; aniline","Industrial"))
      return("Industrial")
    else if(x%in%c("PAH; industrial","PAH"))
      return("PAH")
    else if(x%in%c("personalCare; personalCare; natural; natural; consumerProduct; consumerProduct",
                   "personalCare; natural; consumerProduct","personalCare; natural",
                   "pharmacological; personalCare; industrial; natural; consumerProduct"))
        return("Natural")
    else if(is.na(x)||x=='NA')
        return("Unclassified")
    else return(x))
  return(full.class)

}

#here we join the chemical metadata from the comptox dashboard
chemMeta<-readxl::read_xls(paste0(data.dir,
                                  '/CompToxChemicalsDashboard-Batch-Search_2020-08-27_17_26_49.xls'))%>%
    select(INPUT,DTXSID,PREFERRED_NAME,INCHIKEY,SMILES,MOLECULAR_FORMULA,
           AVERAGE_MASS,PUBCHEM_DATA_SOURCES)%>%
  rename(cas_number='INPUT')%>%
  mutate(chemDescription='WikipediaDataGoesHere')%>%
  right_join(fullMapping)%>%
    left_join(getNewChemicalClass())%>%
    tidyr::replace_na(list(newClass='Unclassified'))%>%
    select(-c(chemical_class,Classification))%>%
    rename(chemical_class='newClass')

chemMeta<-mutate(chemMeta,Chemical_ID=unlist(as.character(zf.cid)))

#These files have the extract information and the mappin gof extracts to chemicals
sampMeta<-readxl::read_xlsx(paste0(data.dir,'/superfund_location_data.xlsx'))
sampChem<-read.csv(paste0(data.dir,'/curatedSampData.csv'))%>%
    select(chem_lims_num,date_sampled,sample_matrix,technology,Sample_ID='zaap_cid',
           zf_lims_id,measurement_value,unit,cas_number)%>%distinct()%>%
  rename(SampleNumber='chem_lims_num')%>%left_join(sampMeta)

#'combineChemicalEndpointData
#'@param bmdfiles a list of files that come from the BMD pipeline
#'@param is_extract - if our chemical id is an extract we use a different metadata
#'@return a data.frame 
##We will release an 'endpoint file for each condition'
combineChemicalEndpointData<-function(bmdfiles,is_extract=FALSE){

  files = lapply(bmdfiles,read.csv)
  
  mid.bmd<-do.call(rbind,files)%>%
#    filter(!is.null(BMD_Analysis_Flag))%>%
    dplyr::select(Chemical_ID,End_Point,Model,BMD10,BMD50,AUC_Norm,BMD10_Flag,BMD50_Flag)#%>%
   # mutate(endPointLink='')
    
  
  if(is_extract)
    full.bmd<-mid.bmd%>%mutate(Sample_ID=as.character(Chemical_ID))%>%
      left_join(sampChem)%>%#%>%mutate(Chemical_ID=as.character(zaap_cid)))%>%
      left_join(endpointDetails)%>%
      select(-c(Chemical_ID,End_Point,cas_number,measurement_value,unit))%>%
    distinct()
  
  else
    full.bmd<-mid.bmd%>%
    mutate(Chemical_ID=as.character(Chemical_ID))%>%
      left_join(chemMeta)%>%
#      rename(Chemical_ID='zf.cid')%>%
      left_join(endpointDetails)%>%
      distinct()%>%select(-c(End_Point))##should we remove endpoint YES

  return(full.bmd)
}


#'combineChemicalFit files 
#'@param bmdfiles the output of files from teh BMD pipeline describing the dose response 
# @return data.frame
##We will release an 'endpoint file for each condition'
combineChemicalFitData<-function(bmdfiles, is_extract=FALSE){

  files = lapply(bmdfiles,read.csv)
  
  full.bmd<-do.call(rbind,files)%>%
    mutate(zf.cid=as.character(Chemical_ID))%>%
    rename(ChemicalId='zf.cid')%>%
     left_join(endpointDetails)%>%
     distinct()%>%
    select(-c(End_Point,Description,ChemicalId))
  if(is_extract)
    full.bmd <- rename(full.bmd,Sample_ID='Chemical_ID')
  return(full.bmd)
}

combineChemicalDoseData<-function(bmdfiles, is_extract=FALSE){

  files = lapply(bmdfiles,read.csv)
  full.bmd<-do.call(rbind,files)%>%
    mutate(zf.cid=as.character(Chemical_ID))%>%
    rename(ChemicalId='zf.cid')%>%
   left_join(endpointDetails)%>%
    dplyr::select(-c(End_Point,Description))%>%
    distinct()%>%select(-ChemicalId)
  if(is_extract)
    full.bmd <- rename(full.bmd,Sample_ID='Chemical_ID')

  return(full.bmd)
}

###ADD IN NEW DATA FILES FROM PARITOSH HERE
#2020-08-06 files here:
bmd.files<-c(paste0(data.dir,'/bmd_vals_morphological_sara.csv'),
             paste0(data.dir,'/bmd_vals_behavioral_sara.csv'),
             paste0(data.dir,'/bmd_vals_NC24_sara.csv')) #added 9/18
curv.files<-c(paste0(data.dir,'/fit_vals_behavioral (1).csv'),
              paste0(data.dir,'/fit_vals_morphological.csv'),
              paste0(data.dir,'/fit_vals_NC24.csv')) #added 9/18
dose.files<-c(paste0(data.dir,'/dose_response_vals_morphological.csv'),
              paste0(data.dir,'/dose_response_vals_behavioral (1).csv'),
              paste0(data.dir,'/dose_response_vals_NC24.csv')) #added 9/18

#additional files processed and sent 9/2/2020
bmd.files <-c(bmd.files,paste0(data.dir,'/bmd_vals_extra_ep_sara.csv'))
curv.files <-c(curv.files,paste0(data.dir,'/fit_vals_extra_ep.csv'))
dose.files <-c(dose.files,paste0(data.dir,'/dose_response_vals_extra_ep.csv'))


#extract data
e.bmd<-paste0(data.dir,'/bmd_vals_extracts_sara.csv')
e.curv<-paste0(data.dir,'/fit_vals_extracts.csv')
e.dose<-paste0(data.dir,'/dose_response_vals_extracts.csv')


##now we check for additional files
parser = argparse::ArgumentParser()#'Command line tool to add extract-dose-response data to SRP Analytics Portal')
parser$add_argument('--bmd',dest='new.bmd',type='character',default='',help='BMD file to add to portal')
parser$add_argument('--doseResponse',dest='new.dose',type='character',default='',help='Dose response curve points to add to portal')
parser$add_argument('--coords',dest='new.curv',type='character',default='',help='New curve fit coordinates to plot in SRP analytics portal ')


args=parser$parse_args()
if(args$new.bmd!="")
    bmd.files<-c(args$new.bmd,bmd.files)
if(args$new.dose!="")
    dose.files<-c(dose.files,args$new.dose)
if(args$new.curv!="")
    curv.files<-c(curv.files,args$new.curv)

bmds<-combineChemicalEndpointData(bmd.files)
curves <-combineChemicalFitData(curv.files)
doseReps <-combineChemicalDoseData(dose.files)

ebmds<-combineChemicalEndpointData(c(e.bmd),TRUE)
ecurves <- combineChemicalFitData(c(e.curv),TRUE)
edrs<- combineChemicalDoseData(c(e.dose),TRUE)

id_mapping <-bmds%>%select(Chemical_ID,cas_number,AVERAGE_MASS)%>%distinct()

sampChem <-sampChem%>%left_join(id_mapping)%>%
  mutate(unit=stringr::str_replace_all(unit,'ug/ml','pg/uL'))%>%
  mutate(unit=stringr::str_replace_all(unit,'pg/ul','pg/uL'))%>%
  subset(unit!='ng PAH/ÂµL DMSO')%>%
  rename(concentration=measurement_value)%>%
  mutate(AVERAGE_MASS=as.numeric(AVERAGE_MASS))%>%rowwise()%>%
  mutate(measurement_value=concentration*1000/AVERAGE_MASS)%>%
  dplyr::select(-c(concentration,AVERAGE_MASS,unit))



##Final output for the platform team is these 4 files
write.csv(bmds,file=paste0('chemSummaryStats.csv'),row.names = FALSE)
write.csv(ebmds,file=paste0('envSampSummaryStats.csv'),row.names=FALSE)

write.csv(curves,file=paste0('chemXYcoords.csv'),row.names = FALSE)
write.csv(ecurves,file=paste0('envSampXYcoords.csv'),row.names = FALSE)

write.csv(doseReps,file=paste0('chemdoseResponseVals.csv'),row.names = FALSE)
write.csv(edrs,file=paste0('envSampdoseResponseVals.csv'),row.names = FALSE)


write.csv(sampChem,file=paste0(out.dir,'chemicalsByExtractSample.csv'),row.names=FALSE)

##TODO: zip them up into package together with a readme/description for datahub 
allfiles<-c('README.md',list.files(path='.')[grep('csv',list.files(path='.'))])
print(allfiles)
tar(paste0(out.dir,'srpAnalyticsCompendium.tar.gz'),files=allfiles,compression='gzip')

