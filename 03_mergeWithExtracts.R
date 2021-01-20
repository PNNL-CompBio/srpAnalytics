##now that we have behavioral data and Serberus platform this
##script is designed to
##process everything into a digestable data

require(dplyr)
require(readxl)
require(argparse)
require(xml2)

##The data release will be comprised of 7 files
#' 1- list of environmental samples, their metadata, and the chemical composition (curated sample data)
#' 2- ZF summary statistics for chemicals
#' 3- ZF points to plot for chemicals
#' 4- ZF curves to plot for chemicals
#' 5- ZF summary statistics for samples
#' 6- ZF points to plot for samples
#' 7- ZF curves to plot for samples


                                        #These pathways refer to absolute pathways in the docker image
##setting these three parameters, can be appended
data.dir<-'/srpAnalytics/data/'

required_sample_columns<-c("SampleNumber","date_sampled","sample_matrix","technology",
    "Sample_ID","zf_lims_id","cas_number","ClientName","SampleName","LocationLat","LocationLon",
    "LocationName","LocationAlternateDescription","AlternateName","Chemical_ID","measurement_value",
    "water_concentration","water_concentration_unit")
#extra columns: "bioassay_sample_from_lims","parent_samples_from_lims","child_samples_from_lims")

required_bmd_columns<-list(bmd=c('Chemical_ID','End_Point','Model','BMD10','BMD50',"Min_Dose","Max_Dose",
                                "AUC_Norm","DataQC_Flag","BMD_Analysis_Flag"),#,"BMD10_Flag","BMD50_Flag"),
                          doseRep=c("Chemical_ID","End_Point","Dose","Response","CI_Lo","CI_Hi"),
                          fitVals=c("Chemical_ID","End_Point","X_vals","Y_vals"))

##output directory is fixed
out.dir<-'/tmp/'

#' Get chemical metadata, which is stored in `data.dir`
#' @param data.dir path to standardized data that enables matching across datasets
#' @return data.frame
getChemMetadata<-function(data.dir){

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

                                        #here we join the chemical metadata from the comptox dashboard
    chemMeta<-readxl::read_xls(paste0(data.dir,
                                      'CompToxChemicalsDashboard-Batch-Search_2020-08-27_17_26_49.xls'))%>%
        select(INPUT,DTXSID,PREFERRED_NAME,INCHIKEY,SMILES,MOLECULAR_FORMULA,
               AVERAGE_MASS,PUBCHEM_DATA_SOURCES)%>%
        rename(cas_number='INPUT')%>%
        #TODO: update this once we get info from Lindsey
        mutate(chemDescription='WikipediaDataGoesHere')%>%
        right_join(fullMapping)%>%
        left_join(getNewChemicalClass(data.dir))%>%
        tidyr::replace_na(list(newClass='Unclassified'))%>%
        select(-c(chemical_class,Classification))%>%
        rename(chemical_class='newClass')

    chemMeta<-mutate(chemMeta,Chemical_ID=unlist(as.character(zf.cid)))

    return(chemMeta)
}

#' getEndpointMetadata
#' @param data.dir directory where standard endpoint data is stored
#' @return data.frame
getEndpointMetadata<-function(data.dir){
                                        #here is our pre-defined dictionary
    endpointDetails<-readxl::read_xlsx(paste0(data.dir,'SuperEndpoint Mapping 2020Sept2.xlsx'),
                                       sheet='Dictionary')%>%
        rename(End_Point='Abbreviation',`End Point Name`='Simple name (<20char)')%>%
        rename(endPointLink='Ontology Link')
    return(endpointDetails)
}

#'getNewChemicalClass
#'This file reads in the files and processes the chemical class names
#'to be friendly for the website
#'@param data.dir
#'@return data.frame
getNewChemicalClass<-function(data.dir){
  pahs<-readxl::read_xlsx(paste0(data.dir,'/PAH_and_1530_SRP_Summary.xlsx'),
                         sheet='Master summary-PAHs')%>%
      dplyr::select(cas_number='casrn')%>%
      mutate(Classification='PAH')

  extras<-data.frame(cas_number=c("3074-03-01","7496-02-08","6373-11-01"),
                     Classification='PAH')

  non.pahs<-readxl::read_xlsx(paste0(data.dir,'/PAH_and_1530_SRP_Summary.xlsx'),
                             sheet='Master Summary-Non-PAHs')%>%
    dplyr::select(cas_number='casrn',Classification='classification')

  full.class<-rbind(pahs,extras,non.pahs)
  full.class$newClass<-sapply(full.class$Classification,function(x)
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

#' buildSampleData - takes the curated information and selects the data we need
#' @param data.dir
#' @return data.frame
buildSampleData<-function(data.dir,chemMeta){
  ##New data provided by michael
  sampChem<-read.csv(paste0(data.dir,'/pnnl_bioassay_sample_query_1-14-2021.csv'))%>%
        subset(SampleNumber!='None')%>%
    mutate(Chemical_ID=as.character(Chemical_ID))%>%
        distinct()

  chemDat<-chemMeta%>%
    select(Chemical_ID,cas_number,AVERAGE_MASS,PREFERRED_NAME,chemical_class)%>%
    distinct()%>%
    mutate(Chemical_ID=as.character(Chemical_ID))

  print("Updating concentration data")
  sampChem <-sampChem%>%
    left_join(chemDat)%>%
 #   mutate(unit<-stringr::str_replace_all(unit,'ug/ml','pg/uL'))%>%
#    mutate(unit<-stringr::str_replace_all(unit,'pg/ul','pg/uL'))%>%
#    subset(unit!='ng PAH/ÂµL DMSO')%>%
    rename(concentration=measurement_value)%>%
    mutate(AVERAGE_MASS=as.numeric(AVERAGE_MASS))%>%
    rowwise()%>%
    mutate(measurement_value=concentration*1000/AVERAGE_MASS)%>%
    dplyr::select(-c(concentration,AVERAGE_MASS))

    return(sampChem)
 }

#'combineChemicalEndpointData
#'@param bmdfiles a list of files that come from the BMD pipeline
#'@param is_extract - if our chemical id is an extract we use a different metadata
#'@return a data.frame
##We will release an 'endpoint file for each condition'
combineChemicalEndpointData<-function(bmdfiles,is_extract=FALSE,sampChem,endpointDetails){

    print(paste('Combining bmd files:',paste(bmdfiles,collapse=',')))
  cols <- required_bmd_columns$bmd
  files <- lapply(bmdfiles,function(x) read.csv(x)%>%dplyr::select(cols))

  mid.bmd<-do.call(rbind,files)%>%
#    filter(!is.null(BMD_Analysis_Flag))%>%
    dplyr::select(cols)

  if(is_extract){
    full.bmd<-mid.bmd%>%
      dplyr::mutate(`Sample_ID`=as.character(Chemical_ID))%>%
      dplyr::select(-Chemical_ID)%>%
      left_join(sampChem,by='Sample_ID')%>%#%>%mutate(Chemical_ID<-as.character(zaap_cid)))%>%
      left_join(endpointDetails)%>%
      dplyr::select(-c('End_Point','cas_number','measurement_value'))%>%
    distinct()
}
  else{
    full.bmd<-mid.bmd%>%
      dplyr::mutate(`Chemical_ID`=as.character(Chemical_ID))%>%
      left_join(sampChem)%>%
#      rename(Chemical_ID<-'zf.cid')%>%
      left_join(endpointDetails)%>%
      distinct()%>%select(-c(End_Point))##should we remove endpoint YES
}
  return(full.bmd)
}


#'combineChemicalFit files
#'@param bmdfiles the output of files from teh BMD pipeline describing the dose response
# @return data.frame
##We will release an 'endpoint file for each condition'
combineChemicalFitData<-function(bmdfiles, is_extract=FALSE, endpointDetails){

    print(paste('Combining fit files:',paste(bmdfiles,collapse=',')))
    files <- lapply(bmdfiles,read.csv)

    full.bmd<-do.call(rbind,files)%>%
        dplyr::select(required_bmd_columns$fitVals)%>%
        mutate(zf.cid=as.character(Chemical_ID))%>%
    rename(ChemicalId='zf.cid')%>%
     left_join(endpointDetails)%>%
     distinct()%>%
    select(-c(End_Point,Description,ChemicalId))
  if(is_extract)
    full.bmd <- rename(full.bmd,Sample_ID='Chemical_ID')
  return(full.bmd)
}

combineChemicalDoseData<-function(bmdfiles, is_extract=FALSE, endpointDetails){
    files <- lapply(bmdfiles,read.csv)
    print(paste('Combining dose response files:',paste(bmdfiles,collapse=',')))
    full.bmd<-do.call(rbind,files)%>%
        dplyr::select(required_bmd_columns$doseRep)%>%
        mutate(zf.cid=as.character(Chemical_ID))%>%
        rename(ChemicalId='zf.cid')%>%
        left_join(endpointDetails)%>%
        dplyr::select(-c(End_Point,Description))%>%
        distinct()%>%select(-ChemicalId)

    if(is_extract)
        full.bmd <- rename(full.bmd,Sample_ID='Chemical_ID')

  return(full.bmd)
}


chem_dirs=c('phase_I_II','zf_morphology')
extract_dirs=c('extracts')

#' main method
#' Parsers arguments
main<-function(){
  ##now we check for additional files
  parser <- ArgumentParser()
  parser$add_argument('-s','--samples',dest='samp_files',default="",
                      help='The subsequent files are samples')
  parser$add_argument('-c','--chemicals',dest='chem_files',default='',help='The subsequent files are chemicals')


  args <- parser$parse_args()
                                        #if we are adding new data, add to additional data in repo

  #files that we're reading in
  chem.files<-unlist(strsplit(args$chem_files,split=','))
  extract.files<-unlist(strsplit(args$samp_files,split=','))

  bmd.files<-c()
  dose.files<-c()
  curv.files<-c()

  e.bmd<-c()
  e.dose<-c()
  e.curve<-c()

    for(chem in chem_dirs){
        path=paste0(data.dir,'/',chem,'/')
        bmd.files<-c(bmd.files,paste0(path,'bmd_vals_all_qc.csv'))
        dose.files<-c(dose.files,paste0(path,'dose_response_vals_all_qc.csv'))
        curv.files<-c(curv.files,paste0(path,'fit_vals_all_qc.csv'))
    }

    for(ext in extract_dirs){
        path=paste0(data.dir,'/',ext,'/')
        e.bmd<-c(e.bmd,paste0(path,'bmd_vals_all_qc.csv'))
        e.dose<-c(e.dose,paste0(path,'dose_response_vals_all_qc.csv'))
        e.curve<-c(e.curve,paste0(path,'fit_vals_all_qc.csv'))
    }


    ##read in files
   if(length(chem.files)==3){
       #path=paste0(data.dir,'/',chem,'/')
       bmd.files<-c(bmd.files,chem.files[1])#paste0(path,'bmd_vals_all_qc.csv'))
       dose.files<-c(dose.files,chem.files[3])#paste0(path)'dose_response_vals_all_qc.csv'
       curv.files<-c(curv.files,chem.files[2])#paste0(path,'fit_vals_all_qc.csv'))
       }
       else{
           print("Not adding any chemical files")
       }
                                        #for(ext in extract.dirs){
    if(length(extract.files)==3){
#        path=paste0(data.dir,'/',ext,'/')
        e.bmd<-c(e.bmd,extract.files[1])#paste0(path,'bmd_vals_all_qc.csv'))
        e.dose<-c(e.dose,extract.files[3])##paste0(path,'dose_response_vals_all_qc.csv'))
        e.curve<-c(e.curve,extract.files[2])#paste0(path,'fit_vals_all_qc.csv'))
    }else{
        print("Not adding any envinromental files")
    }
  chemMeta<-getChemMetadata(data.dir)
  sampChem<-buildSampleData(data.dir,chemMeta)

  print('Processing extract response data')
  endpointDetails<-getEndpointMetadata(data.dir)
  ebmds<-combineChemicalEndpointData(e.bmd,is_extract=TRUE,sampChem,endpointDetails)
  ecurves <- combineChemicalFitData(e.curve,is_extract=TRUE, endpointDetails)
  edrs<- combineChemicalDoseData(e.dose,is_extract=TRUE, endpointDetails)
  print('Processing chemical response data')
  bmds<-combineChemicalEndpointData(bmd.files,FALSE,sampChem,endpointDetails)
  curves <-combineChemicalFitData(curv.files, FALSE, endpointDetails)
  doseReps <-combineChemicalDoseData(dose.files, FALSE, endpointDetails)

  ##there are mismatches, so we should figure out where those exists
  missing<-list(zebrafishNoChem=setdiff(ebmds$Sample_ID,as.character(sampChem$Sample_ID)),
                chemDataNoZebrafish=setdiff(as.character(sampChem$Sample_ID),ebmds$Sample_ID))

  print(missing)
  ##Final output for the platform team is these 4 files
  write.csv(bmds,file=paste0('chemSummaryStats.csv'),row.names = FALSE)
  write.csv(ebmds,file=paste0('envSampSummaryStats.csv'),row.names=FALSE)

  write.csv(curves,file=paste0('chemXYcoords.csv'),row.names = FALSE)
  write.csv(ecurves,file=paste0('envSampXYcoords.csv'),row.names = FALSE)

  write.csv(doseReps,file=paste0('chemdoseResponseVals.csv'),row.names = FALSE)
  write.csv(edrs,file=paste0('envSampdoseResponseVals.csv'),row.names = FALSE)

  write.csv(sampChem,file=paste0('chemicalsByExtractSample.csv'),row.names=FALSE)

  wd <- paste0(getwd(),'/')
    allfiles<-paste0(wd, c('README.md',list.files(path='.')[grep('csv',list.files(path='.'))]))
    print(allfiles)
  print(paste('Now zipping up',length(allfiles),'files'))
  tar(paste0(out.dir,'srpAnalyticsCompendium.tar.gz'),files=allfiles,compression='gzip')
}

main()
