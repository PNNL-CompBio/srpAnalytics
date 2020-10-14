##now that we have behavioral data and Serberus platform this script is designed to
#process everything into a digestable data

require(dplyr)
require(readxl)
##The data release will be comprised of 5 files
#' 1- list of chemicals nad their information
#' 2- list of sites, their metadata, and the chemical composition (curated sample data)
#' 3- ZF summary statistics
#' 4- ZF points to plot
#' 5- ZF curves to plot


#step 1 - get cas ids to class
chemMapping<-read.csv('data/Chemicals.csv')%>%dplyr::select(cas_number,chemical_class,zaap_cid)%>%
  distinct()%>%
  rename(zf.cid='zaap_cid')
newMapping<-readxl::read_xlsx('data/chemicalid_table_2020AUG07.xlsx')%>%
  select(cas_number='casrn',zf.cid='chemical_id')%>%mutate(chemical_class=NA)
fullMapping<-readxl::read_xlsx('data/OSU_BU_Overlap_Inventory_2020MAR30.xlsx')%>%select(casrn,zf.cid)%>%
  distinct()%>%
  tidyr::separate_rows(zf.cid,sep=';\\s+')%>%
  #rename(cas_number='casrn')%>%
  left_join(chemMapping)%>%
  select(-cas_number)%>%
  rename(cas_number='casrn')%>%rbind(newMapping)


endpointDetails<-readxl::read_xlsx('data/SuperEndpoint Mapping 2020Sept2.xlsx', sheet='Dictionary')%>%
  rename(End_Point='Abbreviation',`End Point Name`='Simple name (<20char)')

getNewChemicalClass<-function(){
  pahs=readxl::read_xlsx('data/PAH_and_1530_SRP_Summary.xlsx',sheet='Master summary-PAHs')%>%
    dplyr::select(cas_number='casrn')%>%mutate(Classification='PAH')
  extras<-data.frame(cas_number=c("3074-03-01","7496-02-08","6373-11-01"),Classification='PAH')
  non.pahs=readxl::read_xlsx('data/PAH_and_1530_SRP_Summary.xlsx',sheet='Master Summary-Non-PAHs')%>%
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
    else return(x))
  return(full.class)

}

chemMeta<-readxl::read_xls('data/CompToxChemicalsDashboard-Batch-Search_2020-08-27_17_26_49.xls')%>%
  select(INPUT,DTXSID,PREFERRED_NAME,INCHIKEY,SMILES,MOLECULAR_FORMULA,AVERAGE_MASS,PUBCHEM_DATA_SOURCES)%>%
  rename(cas_number='INPUT')%>%
  right_join(fullMapping)%>%
  left_join(getNewChemicalClass())%>%
  select(-c(chemical_class,Classification))%>%
  rename(chemical_class='newClass')

chemMeta$zf.cid<-unlist(as.character(chemMeta$zf.cid))

#full characterization of sample data with ids
sampMeta<-readxl::read_xlsx('data/superfund_location_data.xlsx')
sampChem<-read.csv('data/curatedSampData.csv')%>%
  select(chem_lims_num,date_sampled,sample_matrix,technology,zaap_cid,extract_exp,zf_lims_id,measurement_value,unit,cas_number)%>%distinct()%>%
  rename(SampleNumber='chem_lims_num')%>%left_join(sampMeta)

##We will release an 'endpoint file for each condition'
combineChemicalEndpointData<-function(bmdfiles){

  files = lapply(bmdfiles,read.csv)
  full.bmd<-do.call(rbind,files)%>%
#    filter(!is.null(BMD_Analysis_Flag))%>%
    dplyr::select(Chemical_ID,End_Point,Model,BMD10,BMD50,AUC_Norm,BMD10_Flag,BMD50_Flag)%>%
    mutate(zf.cid=as.character(Chemical_ID))%>%
    left_join(chemMeta)%>%
    rename(ChemicalId='zf.cid')%>%
    left_join(endpointDetails)

  return(full.bmd)
}


##We will release an 'endpoint file for each condition'
combineChemicalFitData<-function(bmdfiles){

  files = lapply(bmdfiles,read.csv)
  full.bmd<-do.call(rbind,files)%>%
    #filter(!is.null(BMD_Analysis_Flag))%>%
   # dplyr::select(Chemical_ID,End_Point,Model,BMD10,BMDL,BMD50,AUC,Min_Dose,Max_Dose)%>%
    mutate(zf.cid=as.character(Chemical_ID))%>%
 #   left_join(chemMeta,by='zf.cid')%>%
    rename(ChemicalId='zf.cid')%>%
    left_join(endpointDetails)

  return(full.bmd)
}

combineChemicalDoseData<-function(bmdfiles){

  files = lapply(bmdfiles,read.csv)
  full.bmd<-do.call(rbind,files)%>%
    #filter(!is.null(BMD_Analysis_Flag))%>%
    # dplyr::select(Chemical_ID,End_Point,Model,BMD10,BMDL,BMD50,AUC,Min_Dose,Max_Dose)%>%
    mutate(zf.cid=as.character(Chemical_ID))%>%
#    left_join(chemMeta,by='zf.cid%>%')%>%
    rename(ChemicalId='zf.cid')%>%
    left_join(endpointDetails)

  return(full.bmd)
}

###ADD IN NEW DATA FILES FROM PARITOSH HERE
#2020-08-06 files here:
bmd.files<-c('data/bmd_vals_morphological_sara.csv',
             'data/bmd_vals_behavioral_sara.csv',
             'data/bmd_vals_NC24_sara.csv') #added 9/18
curv.files<-c('data/fit_vals_behavioral (1).csv',
              'data/fit_vals_morphological.csv',
              'data/fit_vals_NC24.csv') #added 9/18
dose.files<-c('data/dose_response_vals_morphological.csv',
              'data/dose_response_vals_behavioral (1).csv',
              'data/dose_response_vals_NC24.csv') #added 9/18

#additional files processed and sent 9/2/2020
bmd.files <-c(bmd.files,'data/bmd_vals_extra_ep_sara.csv')
curv.files <-c(curv.files,'data/fit_vals_extra_ep.csv')
dose.files <-c(dose.files,'data/dose_response_vals_extra_ep.csv')

bmds<-combineChemicalEndpointData(bmd.files)
curves <-combineChemicalFitData(curv.files)
doseReps <-combineChemicalDoseData(dose.files)

id_mapping <-bmds%>%select(Chemical_ID,cas_number,AVERAGE_MASS)%>%distinct()
sampChem <-sampChem%>%left_join(id_mapping)%>%
  mutate(unit=stringr::str_replace_all(unit,'ug/ml','pg/uL'))%>%
  mutate(unit=stringr::str_replace_all(unit,'pg/ul','pg/uL'))%>%
  subset(unit!='ng PAH/ÂµL DMSO')%>%
  rename(concentration=measurement_value)%>%
  mutate(AVERAGE_MASS=as.numeric(AVERAGE_MASS))%>%rowwise()%>%
  mutate(measurement_value=concentration*1000/AVERAGE_MASS)%>%
  dplyr::select(-c(concentration,AVERAGE_MASS,unit))


write.csv(bmds,file='chemSummaryStats.csv',row.names = FALSE)
write.csv(curves,file='chemXYcoords.csv',row.names = FALSE)
write.csv(doseReps,file='chemdoseResponseVals.csv',row.names = FALSE)
write.csv(sampChem,file='chemicalsByExtractSample.csv',row.names=FALSE)

#chemsOnly<-bmds%>%dplyr::select(Chemical_ID,cas_number,PREFERRED_NAME,chemical_class)%>%distinct()
#write.csv(chemsOnly,'chemsWithEndpoints.csv',row.names=F)

##do we have curves for the extract data? if so we need to merge that as well
