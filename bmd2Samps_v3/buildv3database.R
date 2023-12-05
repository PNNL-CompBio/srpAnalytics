##V2 database will be comprised of more files to distinctly track
## gene, chemical, and sample metadata

require(dplyr)
require(readxl)
require(argparse)
require(xml2)

##The data release will be comprised of 9 files (note change from v1!)
#' 1- list of environmental samples and the chemical composition (curated sample data)
#' 2- ZF summary statistics for chemicals (and no chemical metadata)
#' 3- ZF points to plot for chemicals
#' 4- ZF curves to plot for chemicals
#' 5- ZF summary statistics for samples (and no sample metadata)
#' 6- ZF points to plot for samples
#' 7- ZF curves to plot for samples
#' 8- new for v2: chemical metdata
#' 9- new for v2: sample metadata

#These pathways refer to absolute pathways in the docker image
##setting these three parameters, can be appended
data.dir<-'/bmd2Samps/data/'

#data.dir='./data/'
##output directory is fixed
out.dir<-'/tmp/'
#out.dir<-'./'



##################################
#Master ID tables
#The database requires Sample_ID and  Chemical_ID be unique. They are in some files but not others
#so i automatically update these tables below
##################################

chemIdMasterTable<-function(casIds){
  map <- read.csv2(paste0(data.dir,'chemicalIdMapping.csv'),header=TRUE,sep=',',
                   fileEncoding="UTF-8-BOM")%>%
    dplyr::select(cas_number,zf.cid,Chemical_ID,chemical_class)%>%
      distinct()%>%
      subset(!is.na(cas_number))

  ###there are some acceptable duplicates:
  ## 3144, 3143 - same cas, different BMD
  ## 1161, 1171
  ## 441,447
  missing <- setdiff(casIds,map$cas_number)
  if(length(missing)>0){
    message("Missing ",length(missing),' chemical IDs, adding them now')
    maxId <- max(map$Chemical_ID)+1
    newMax <- maxId + length(missing)-1
    newMap <- rbind(map, data.frame(cas_number=missing,zf.cid=rep('',length(missing)),
                                    Chemical_ID=seq(maxId,newMax),
                                    chemical_class=rep('',length(missing))))
    write.csv(newMap,paste0(data.dir,'chemicalIdMapping.csv'),row.names = F)
    map <- newMap
  }
  return(map)

}

sampIdMasterTable<-function(existingSampNumbers){

  map <- read.csv(paste0(data.dir,'sampleIdMapping.csv'),fileEncoding="UTF-8-BOM")%>%
    select(Sample_ID,SampleNumber)%>%
    distinct()

  missing <- setdiff(existingSampNumbers,map$SampleNumber)
  if(length(missing)>0){
    message("Missing ",length(missing),' sample IDs, adding them now')
    maxId <- max(as.numeric(map$Sample_ID),na.rm=T)+1 ##some are going to be lost, that's ok.
    newMax <- maxId + length(missing)-1
    newMap <- rbind(map, data.frame(Sample_ID=seq(maxId,newMax),
                                    SampleNumber=missing))
    write.csv(newMap,paste0(data.dir,'sampleIdMapping.csv'),row.names = F)
    map <- newMap
  }
  return(map)
}


###################################
# Duplicate removal
# there are two causes of duplicates in the data
# First, there are samples that are evaluated from multiple files, so we need to remove
# those that are duplicated
# Second, there are multiple chemical ids mapping to a single CAS id. In this case,
# We need to select a single chemical ids
###################################

#' Remove duplicates from anything with a Chemical_ID in the field
#' @param tab Table to be de-duplicated
#' @param chemIds Mapping of chemical_IDs to CAS ids
#' @param cols List of columsn that make the data distinct
#' @return data.frame of de-duplicated table
removeChemIdDuplicates<-function(tab, chemIds,cols=c('AUC_Norm','End_Point_Name','Model')){
  
  #first get all duplicaetd chemical ids
  dupes<-chemIds|>group_by(cas_number)|>summarize(nc=n_distinct(Chemical_ID))|>subset(nc>1)
  chem_ids=subset(chemIds,cas_number%in%dupes$cas_number)|>
    dplyr::select(cas_number,Chemical_ID)|>
    mutate(inDataset=Chemical_ID%in%tab$Chemical_ID)|>
    subset(inDataset)
  
  ##now which ones are actually in the data?

  #create a distinct table
  dist.tab<-tab|>
    left_join(chemIds)|>
    dplyr::select(c('cas_number','Chemical_ID',cols))
  
#  dupe.counts<- dist.tab|>group_by(c(cols,'cas_number'))|>
#    summarize(chemCounts=n_distinct(Chemical_ID))

  dupe.counts<- dist.tab|>
    group_by(cas_number,End_Point_Name,Model)|>
    summarize(nIds=n_distinct(Chemical_ID))|>
    subset(nIds>1)
  
  dupe.cas<-unique(dupe.counts$cas_number)
  
}
###################################
# Metadata collection
# these functions import metadata
###################################

#' Get chemical metadata, which is stored in `data.dir`
#' @param data.dir path to standardized data that enables matching across datasets
#' @param comptoxfiles Files downloaded from the EPA website
#' @return data.frame
getChemMetadata<-function(data.dir,
                          comptoxfiles=c('CompToxChemicalsDashboard-Batch-Search_2021-11-05_17_57_46.xlsx',
                                         'CCD-Batch-Search_2022-01-26_10_28_30.xlsx')){    ##This mapping file consumes data from OSU to match identifiers to CAS

     ##we have curated descriptions for each chemical
   curatedDesc <- readxl::read_xlsx(paste0(data.dir,'ChemicalDescriptions.xlsx'))%>%
        select(cas_number='CASRN',chemDescription='USE CATEGORY/DESCRIPTION')

   ##we have chemical source and class information
   chemicalClasses <- masvChemClass(data.dir)
    ##this file comes from COMPTOX
    ##here we join the chemical metadata from the comptox dashboard
   

    chemMeta<-do.call(rbind,lapply(comptoxfiles,function(x) readxl::read_xlsx(paste0(data.dir,x),
                                                                              sheet='Main Data')%>%
                                      select(required_comptox_columns)))%>%
        rename(cas_number='INPUT')%>%
        left_join(chemicalClasses)%>% ##add in chemical class
        left_join(curatedDesc)%>%
        subset(!is.na(cas_number))%>%
        subset(!cas_number%in%c('NA','N/A'))
        #rename(chemical_class='newClass')
      ##we have a mapping of cas to chemical_ids - add in extra if not there
    
    ##manual list of chemical ids - some from OSU, some I created
    chemIds <- chemIdMasterTable(chemMeta$cas_number)%>%
      select(cas_number,Chemical_ID)%>%distinct()
    
    
    chemMeta <- chemMeta%>%
      left_join(chemIds)%>%
      full_join(chemicalClasses)%>% ## add in description data
      tidyr::replace_na(list(newClass='Unclassified'))%>% ##allow for unclassified
      tidyr::replace_na(list(PREFERRED_NAME='Chemical name unknown'))%>% ##should we call this something else?
      mutate(PREFERRED_NAME=stringr::str_replace_all(PREFERRED_NAME,'^-$',"Chemical name unknown"))%>%
      select(-c(ParameterName))
    
    nocas=grep("NOCAS",chemMeta$cas_number)
    if(length(nocas)>0){
        message(paste0('removing ',length(nocas),' chems with no cas'))
        chemMeta<-chemMeta[-nocas,]
    }
    return(chemMeta)
}

#' getEndpointMetadata
#' @param data.dir directory where standard endpoint data is stored
#' @return data.frame
getEndpointMetadata<-function(data.dir){
                                        #here is our pre-defined dictionary
    endpointDetails<-readxl::read_xlsx(paste0(data.dir,'SuperEndpoint Mapping 2021NOV04.xlsx'),
                                       sheet='Dictionary')%>%
        #subset(`Portal Display`=='Display')%>%
        rename(End_Point='Abbreviation',`End_Point_Name`='Simple name (<20char)')%>%
        rename(endPointLink='Ontology Link')%>%
        #select(-`Portal Display`)%>%
      mutate(End_Point=stringr::str_trim(End_Point))%>%
    ##add in endpoint detail for no data
      rbind(list(IncludeInPortal='No',End_Point='NoData',`End_Point_Name`=NA,Description='No data',endPointLink=''))
    

    return(endpointDetails)
}

#'getNewChemicalClass
#'This file reads in the files and processes the chemical class names
#'to be friendly for the website
#' @depracated As we move to the new MASV classes
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
    if(x%in%c('industrial',"industrial; aniline","Industrial",
    "industrial; consumerProduct; phenol","industrial; consumerProduct; aniline","industrial; phenol" ))
      return("Industrial")
    else if(x%in%c("PAH; industrial","PAH"))
      return("PAH")
    else if(x%in%c("personalCare; personalCare; natural; natural; consumerProduct; consumerProduct",
                   "personalCare; natural; consumerProduct","personalCare; natural",
                   "pharmacological; personalCare; industrial; natural; consumerProduct"))
        return("Natural")
    else if(x=='pestFungicide')
      return('Fungicide')
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
    sampChem <- read.csv(paste0(data.dir,'/fses/fses_data_for_pnnl_4-27-2021.csv'),fileEncoding="UTF-8-BOM")%>%
                                        #  sampChem<-read.csv(paste0(data.dir,'/pnnl_bioassay_sample_query_1-14-2021.csv'))%>%
        dplyr::select(required_sample_columns)%>%
        subset(SampleNumber!='None')%>%
        subset(cas_number!='NULL')%>%
        mutate(water_concentration_molar=stringr::str_replace_all(water_concentration_molar,'BLOD|NULL|nc:BDL',"0"))%>%
        mutate(measurement_value_molar=stringr::str_replace_all(measurement_value_molar,'BLOD|NULL|BDL',"0"))%>%
        mutate(water_concentration=stringr::str_replace_all(water_concentration,'BLOD|NULL|BDL',"0"))%>%
                                        # subset(water_concentration_molar!='0.0')%>%
        subset(!measurement_value_molar%in%c('0'))%>%
        subset(!measurement_value%in%c("0","NULL",""))#%>%
#        select(-c(Sample_ID))#,Chemical_ID)) ##These two are added in the 4/27 version of the file

    ##data added 1/19/2022
    newSamp <- readxl::read_xlsx(paste0(data.dir,'/fses/FSES_indoor_outdoor_study.xlsx'))%>%
        dplyr::select(required_sample_columns)%>%
      mutate(LocationLon=LocationLon*-1)

    chemDat<-chemMeta%>%
        select(Chemical_ID,cas_number,AVERAGE_MASS)%>%#,PREFERRED_NAME,chemDescription)%>%
        distinct()

    finalSampChem <-sampChem%>%
        rbind(newSamp)

    ##This mapipng file maps tanguay lab identifiers to those in the anderson lab
    ids<-sampIdMasterTable(finalSampChem$SampleNumber)

    finalSampChem <- finalSampChem %>%
      left_join(chemDat,by='cas_number')%>%
      left_join(ids,by='SampleNumber')%>%
      distinct()
    ##now we have to fix duplicate sample names

    #get duplicates, assign
    all.samp.names<-finalSampChem%>%
      select(Sample_ID,SampleName)%>%
      distinct()%>%
      mutate(isDupe=duplicated(SampleName))

    #filter only duplicates
    dupe.samp.names<-all.samp.names%>%
      subset(isDupe)%>%
      arrange(SampleName)
    
    #get number
    num.dupes<-dupe.samp.names%>%
      group_by(SampleName)%>%
      summarize(nid=n_distinct(Sample_ID))%>%
      left_join(dupe.samp.names)
    
    #paste sequence to end of each name
    new.names<-num.dupes%>%
      select(SampleName,nid)%>%
      distinct()%>%
      rowwise()%>%
      mutate(newName=paste(SampleName,seq(nid),collapse=':'))%>%
      tidyr::separate_rows(newName,sep=':')%>%
      select(oldSN='SampleName',newName)%>%
      arrange(oldSN)
    
    #combine in new data frame
    full.rep=data.frame(dupe.samp.names,new.names)%>%
      select(Sample_ID,SampleName='newName')

    new.samp.names<-all.samp.names%>%
      subset(!isDupe)%>%
      select(-isDupe)%>%
      rbind(full.rep)

    finalSampChem<-finalSampChem%>%
      select(-SampleName)%>%#remove duplicates
      left_join(new.samp.names)

    message("Updating concentration data")
                                        #TODO: separate those without molar concentrations, recompute
    blanks<-finalSampChem%>%subset(measurement_value_molar=='')
    nonblanks<-finalSampChem%>%subset(measurement_value_molar!='')
    fixed<-blanks%>%
        mutate(concentration=as.numeric(measurement_value))%>%
        mutate(AVERAGE_MASS=as.numeric(AVERAGE_MASS))%>%
        rowwise()%>%
        mutate(measurement_value_molar=concentration*1000/AVERAGE_MASS)%>%
        select(-concentration)

    finalSampChem<-rbind(nonblanks,fixed)%>%
      subset(!is.na(cas_number))%>%
        select(-AVERAGE_MASS)%>%
        distinct()

    ##now we have one more rename of samples and metadata
    sampleNameRemap<-readxl::read_xlsx(paste0(data.dir,'/envSampCleanMapping.xlsx'))%>%
      dplyr::select(Sample_ID,date_sampled,sample_matrix,technology,#projectName='ProjectName',SampleName='NewSampleName',
                    #LocationName='NewLocationName')%>%
                    ProjectName,NewSampleName,NewLocationName)%>%
      distinct()

    finalSampChem<-finalSampChem%>%
    #  select(-c(sample_matrix,technology,SampleName,LocationName,date_sampled))%>% ### remove old names
        #select(-c(sample_matrix,technology,date_sampled))%>% ### remove old names
      left_join(sampleNameRemap,by='Sample_ID')%>% ##add in updates
      distinct()
    
    ##older files will have a missing projectName and require manual mapping of some terms
    nas <- which(is.na(finalSampChem$projectName))  ##these are the older samples
    finalSampChem$projectName[nas]<-finalSampChem$ProjectName[nas]
    finalSampChem$LocationName[nas]<-finalSampChem$NewLocationName[nas]
    finalSampChem$SampleName[nas]<-finalSampChem$NewSampleName[nas]
    finalSampChem$date_sampled.x[nas]<-finalSampChem$date_sampled.y[nas]
    finalSampChem$sample_matrix.x[nas]<-finalSampChem$sample_matrix.y[nas]
    finalSampChem$technology.x[nas]<-finalSampChem$technology.y[nas]
    
    finalSampChem<-select(finalSampChem,-c(ProjectName,NewSampleName,NewLocationName,date_sampled.y,sample_matrix.y,technology.y))%>%
      rename(sample_matrix='sample_matrix.x',date_sampled='date_sampled.x',technology='technology.x')%>%
      distinct()%>%
      subset(cas_number!='N/A')
    
    
    ##last thing
    
    return(finalSampChem)

}
#'combineChemicalEndpointData produces the summary statistics from the BMD analysis
#'@param bmdfiles a list of files that come from the BMD pipeline
#'@param is_extract - if our chemical id is an extract we use a different metadata
#'@return a data.frame
##We will release an 'endpoint file for each condition'
# combineChemicalEndpointData<-function(bmdfiles,is_extract=FALSE,sampChem,endpointDetails){
# 
#   ##read in the BMD files formatted by the zf module
#   print(paste('Combining bmd files:',paste(bmdfiles,collapse=',')))
#   cols <- required_bmd_columns$bmd
#   files <- lapply(bmdfiles,function(x) read.csv(x)%>%dplyr::select(cols))
# 
#   mid.bmd<-do.call(rbind,files)%>%
#       dplyr::select(cols)
# 
#   ##some of the chemicals have BMDs that were computed twice and i'm not sure which
#   ##are which anymore. as such, i incorporate the files in chronological order
#   ##and remove the second values
#   dupes<-which(mid.bmd%>%select(Chemical_ID,End_Point)%>%duplicated())
#   if(length(dupes)>0){
#     mid.bmd<-mid.bmd[-dupes,]
#   }
#   
#   if(is_extract){
#     sdSamp<-sampChem%>%
#       tidyr::separate('Sample_ID',into=c('tmpId','sub'),sep='-',remove=FALSE)%>%
#       select(-sub)
# 
#     full.bmd<-mid.bmd%>%
#       dplyr::mutate(tmpId=as.character(Chemical_ID))%>%
#       dplyr::select(-Chemical_ID)%>%
#       full_join(sdSamp,by='tmpId')#%>%#%>%mutate(Chemical_ID<-as.character(zaap_cid)))%>%
# 
#     #fix up sample ids
#     nas<-which(is.na(full.bmd$Sample_ID))
#     full.bmd$Sample_ID[nas]<-full.bmd$tmpId[nas]
# 
#     #now fix up sample names
#     new.nas<-which(is.na(full.bmd$SampleName))
#     full.bmd$SampleName[new.nas]<-paste('Sample',full.bmd$Sample_ID[new.nas])
# 
#     full.bmd<-full.bmd%>%
#       tidyr::replace_na(list(End_Point='NoData'))%>%
#       right_join(endpointDetails)%>%
#       dplyr::select(-c('End_Point','tmpId'))%>%
#       subset(!is.na(Sample_ID))%>%
#     distinct()%>%
#       tidyr::replace_na(list(LocationName='None'))
#   }
#   else{
#     full.bmd<-mid.bmd%>%
#       #dplyr::mutate(`Chemical_ID`=as.character(Chemical_ID))%>%
#       full_join(sampChem)%>%
#       tidyr::replace_na(list(End_Point='NoData'))%>%
# #      rename(Chemical_ID<-'zf.cid')%>%
#       right_join(endpointDetails)%>%
#       subset(!is.na(cas_number))%>%
#       distinct()%>%select(-c('End_Point'))%>%
#       tidyr::replace_na(list(chemical_class='Unclassified'))##should we remove endpoint YES
#   }
# 
#   ##now we fix QC values
#   full.bmd <- full.bmd%>%
#     mutate(DataQC_Flag=ifelse(qc_num%in%c(0,1),'Poor',ifelse(qc_num%in%c(4,5),'Moderate','Good')))%>%
#     rowwise()%>%
#     mutate(Model=stringr::str_replace_all(Model,"NULL","None"))%>%
#       select(-c(qc_num,BMD_Analysis_Flag))
#   
# 
# 
#   return(full.bmd)
# }


#'combineChemicalEndpointDataV2 produces the summary statistics from the BMD analysis
#'@param bmdfiles a list of files that come from the BMD pipeline
#'@param is_extract - if our chemical id is an extract we use a different metadata
#'@return a data.frame
##We will release an 'endpoint file for each condition'
combineV2ChemicalEndpointData<-function(bmdfiles,is_extract=FALSE,sampChem,endpointDetails){
  
  ##read in the BMD files formatted by the zf module
  print(paste('Combining bmd files:',paste(bmdfiles,collapse=',')))
  cols <- required_bmd_columns$bmd
  files <- lapply(bmdfiles,function(x) read.csv(x)%>%dplyr::select(cols))
  
  mid.bmd<-do.call(rbind,files)
  
  ##some of the chemicals have BMDs that were computed twice and i'm not sure which
  ##are which anymore. as such, i incorporate the files in chronological order
  ##and remove the second values
  dupes<-which(mid.bmd%>%select(Chemical_ID,End_Point)%>%duplicated())
  
  if(length(dupes)>0){
    mid.bmd<-mid.bmd[-dupes,]
  }
  
  if(is_extract){
    sdSamp<-sampChem%>%
      tidyr::separate('Sample_ID',into=c('tmpId','sub'),sep='-',remove=FALSE)%>%
      #select(-sub)
      select(tmpId,Sample_ID)%>%distinct()
    
    full.bmd<-mid.bmd%>%
      dplyr::mutate(tmpId=as.character(Chemical_ID))%>%
      dplyr::select(-Chemical_ID)%>%
      full_join(sdSamp,by='tmpId')#%>%#%>%mutate(Chemical_ID<-as.character(zaap_cid)))%>%
    
    #fix up sample ids
    nas<-which(is.na(full.bmd$Sample_ID))
    full.bmd$Sample_ID[nas]<-full.bmd$tmpId[nas]
    
    #now fix up sample names
    new.nas<-which(is.na(full.bmd$SampleName))
    if(length(new.nas)>0)
      full.bmd$SampleName[new.nas]<-paste('Sample',full.bmd$Sample_ID[new.nas])
    
    full.bmd<-full.bmd%>%
      tidyr::replace_na(list(End_Point='NoData'))%>%
      right_join(endpointDetails)%>%
      dplyr::select(-c('End_Point','tmpId'))%>%
      subset(!is.na(Sample_ID))%>%
      distinct()%>%
      tidyr::replace_na(list(LocationName='None'))
  }
  else{
    full.bmd<-mid.bmd%>%
    #V2 update:  full_join(sampChem)%>%
      tidyr::replace_na(list(End_Point='NoData'))%>%
      right_join(endpointDetails)%>%
      #V2 udpate: subset(!is.na(cas_number))%>%
      distinct()%>%
      select(-c('End_Point'))%>%
      tidyr::replace_na(list(chemical_class='Unclassified'))##should we remove endpoint YES
  }
  
  ##now we fix QC values
  full.bmd <- full.bmd%>%
    rename(qc_num='DataQC_Flag')%>%
    mutate(DataQC_Flag=ifelse(qc_num%in%c(0,1),'Poor',ifelse(qc_num%in%c(4,5),'Moderate','Good')))%>%
    rowwise()%>%
    mutate(Model=stringr::str_replace_all(Model,"NULL","None"))%>%
    select(-c(qc_num))
  
  
  
  return(full.bmd)
}


#'combineChemicalFit combines the xy coordinates of the best model fit
#'@param bmdfiles the output of files from the BMD pipeline describing xy coordinates
# @return data.frame
##We will release an 'endpoint file for each condition'
combineChemicalFitData<-function(bmdfiles, is_extract=FALSE, sampChem, endpointDetails){

    print(paste('Combining fit files:',paste(bmdfiles,collapse=',')))
    cols <- required_bmd_columns$fitVals
    files <- lapply(bmdfiles,function (x) read.csv(x)%>%dplyr::select(cols))

    ##get chemicals and EPs in each file, so we can only keep the newest ones
    chemEps<-lapply(files,function(x){
      res<-x%>%select(Chemical_ID,End_Point)%>%
        distinct()%>%
        mutate(combined=paste(Chemical_ID,End_Point))
      return(res$combined)})


    ##now iterate and remove those from earlier files (listed later)
    ##to avoid duplicaates
    newChemEps=chemEps
    ##for loop is messy but best i can do with set diff function

    if(length(chemEps)>1){
        for(i in 2:length(chemEps)){
            orig=c()
            for(j in 1:(i-1))
                orig=union(orig,chemEps[[j]])
            newChemEps[[i]]<-setdiff(newChemEps[[i]],orig)
        }
    }
#    print(newChemEps)
    fixed.files<-lapply(1:length(files),function(i){
        #print(newChemEps[[i]])
        f<-files[[i]]%>%
            mutate(combined=paste(Chemical_ID,End_Point)) #get common index
      f%>%
       subset(combined%in%newChemEps[[i]])%>% ##filter out those that we want
        dplyr::select(cols)%>%
        mutate(zf.cid=as.character(Chemical_ID))%>% ##why did i do this?
        rename(ChemicalId='zf.cid')%>%
        subset(X_vals!="NULL")%>%
        mutate(X_vals=as.numeric(X_vals))%>%
        mutate(Y_vals=as.numeric(Y_vals))

    })

    mid.bmd<-do.call(rbind,fixed.files)

    print(head(mid.bmd))
    full.bmd<-mid.bmd%>%
        right_join(endpointDetails)%>%
        distinct()%>%
        select(-c(End_Point,Description,ChemicalId))
    
    ###TODO: update for reduce data...
    if(is_extract){

       ##make sure we get the information for the sample id
      sdSamp<-sampChem%>%
          tidyr::separate('Sample_ID',into=c('tmpId','sub'),sep='-',remove=FALSE)%>%
      select(Sample_ID,tmpId)%>%
      distinct()

      ##now we join the data withi the updated sample information
    full.bmd<-full.bmd%>%
      dplyr::mutate(tmpId=as.character(Chemical_ID))%>%
        dplyr::select(-Chemical_ID)%>%
        full_join(sdSamp,by='tmpId')%>%
        select(-tmpId)#%>%#%>%mutate(Chemical_ID<-as.character(zaap_cid)))%>%

    #full.bmd <- rename(full.bmd,Sample_ID='Chemical_ID')
    }else{
        full.bmd<-full.bmd%>%
            subset(Chemical_ID%in%sampChem$Chemical_ID)
    }
    return(full.bmd)
}


#' masvChemClass
#' Reads in full MASV class annotations and assigns values to chemicals
#' @return data.frame
masvChemClass<-function(data.dir){
  library(readxl)
  library(tidyr)
  library(dplyr)
  classes<-readxl::read_xlsx(paste0(data.dir,'MASV Classifications 2021.xlsx'),
                             sheet = 'MASV15 Classifications')
  source=c("pharmacological","personalCare","industrial","pulpAndPaper","pestProduct","natural",'pestRodenticide',
           "consumerProduct","pestHerbicide","pestInsecticide",'pestFungicide','pestGeneral','flameRetardant')

  cc=c('PAH','OPAH','PCB','PBB','PBDE','deuterated','dioxinsAndFurans','haloEthers','OPFR','phenol','aniline','uncategorized')


  ##first let's handle  the sources
  chemSource<-classes%>%
    dplyr::select(-cc)%>%
    pivot_longer(cols=source,names_to='chem_source',values_to='posNeg')%>%
  subset(posNeg!='NULL')%>%
  dplyr::select(-posNeg)%>%
    mutate(chem_source=gsub('pest.*','pesticide',chem_source))

  chemSource<-aggregate(chem_source~CASNumber+ParameterName,chemSource,function(x) paste0(x,collapse=';'))


  chemClass <- classes%>%pivot_longer(cols=cc,names_to='chemical_class',values_to='posNeg')%>%
    dplyr::select(-source)%>%
    subset(posNeg!='NULL')%>%
   dplyr::select(-posNeg)

  chemClass<-aggregate(chemical_class~CASNumber+ParameterName,chemClass,function(x) paste0(x,collapse=';'))

  chemComb <- chemSource%>%full_join(chemClass)%>%
    replace_na(list(chemical_class='Unclassified'))%>%
    rename(cas_number='CASNumber')

  write.csv(chemComb,'MASV_classAndSource.csv')
  return(chemComb)

}


#'combineChemicalDoseData combines the dose response points and joins them with metadata
#'@param bmdfiles the output of files in the BMD pipeline describing the dose data
#'@param is_extract Set to true when the samples represent extracts rather than chemicals
#'@param endpointDetails metadata for endpoints
#'@return data.frame
combineChemicalDoseData<-function(bmdfiles, is_extract=FALSE, sampChem,endpointDetails){

    cols<-required_bmd_columns$doseRep
    files <- lapply(bmdfiles,function(x) read.csv(x)%>%select(cols))

    print(paste('Combining dose response files:',paste(bmdfiles,collapse=',')))

    ##get chemicals and EPs in each file, so we can only keep the newest ones
    chemEps<-lapply(files,function(x){
      res<-x%>%select(Chemical_ID,End_Point)%>%
        distinct()%>%
        mutate(combined=paste(Chemical_ID,End_Point))
      return(res$combined)})

    ##now iterate and remove those from earlier files (listed later)
    ##to avoid duplicaates
    newChemEps=chemEps
    
    if(length(chemEps)>1){
      ##for loop is messy but best i can do with set diff function
      for(i in 2:length(chemEps)){
        orig=c()
        for(j in 1:(i-1))
          orig=union(orig,chemEps[[j]])
        newChemEps[[i]]<-setdiff(newChemEps[[i]],orig)
      }
    }

    fixed.files<-lapply(1:length(files),function(i){
      files[[i]]%>%
        mutate(combined=paste(Chemical_ID,End_Point))%>% #get common index
        subset(combined%in%newChemEps[[i]])%>% ##filter out those that we want
        dplyr::select(required_bmd_columns$doseRep)%>%
         mutate(zf.cid=as.character(Chemical_ID))%>% ##why did i do this?
       rename(ChemicalId='zf.cid')
    })

    mid.bmd<-do.call(rbind,fixed.files)

    full.bmd<-mid.bmd%>%
        right_join(endpointDetails)%>%
        dplyr::select(-c(End_Point,Description))%>%
        distinct()%>%select(-ChemicalId)

    if(is_extract){
      sdSamp<-sampChem%>%tidyr::separate('Sample_ID',into=c('tmpId','sub'),sep='-',remove=FALSE)%>%
        select(Sample_ID,tmpId)%>%
        distinct()

      full.bmd<-full.bmd%>%
        dplyr::mutate(tmpId=as.character(Chemical_ID))%>%
        dplyr::select(-Chemical_ID)%>%
          full_join(sdSamp,by='tmpId')%>%
                                        select(-tmpId)#%>%#%>%mutate(Chemical_ID<-as.character(zaap_cid)))%>%

    }#else{
    #    full.bmd<-full.bmd%>%
    #    subset(Chemical_ID%in%sampChem$Chemical_ID)
       # full.bmd <- rename(full.bmd,Sample_ID='Chemical_ID')
  # }
  return(unique(full.bmd))
}


#'
#' buildDB
#' Master function that assembles all the files needed for hte database
#' @param chem.files new chemical bmd files
#' @param extract.files new extract bmd files
buildDB<-function(chem.files=c(),extract.files=c()){
  
  ##these are the subdirectories where we store the original data (Small enough for github)
  chem_dirs=c('phase_I_II_III')#,'zf_morphology')
  extract_dirs=c('extracts')
  
  
  ##here are the files we need in the end - bmds, dose response, and curv fit files for chemicals and extracts
  bmd.files<-c()
  dose.files<-c()
  curv.files<-c()
  
  e.bmd<-c()
  e.dose<-c()
  e.curve<-c()
  
  ##first we collect the files that are output from the bmd code
  for(chem in chem_dirs){
    path=paste0(data.dir,'/',chem,'/')
    bmd.files<-c(bmd.files,paste0(path,c('bmd_vals_all_qc.csv',
                                         'bmd_vals_2021_05_18_all_phase_III_morpho.csv',
                                         'bmd_vals_2021_04_26.csv',
                                         'bmd_vals_2021_01_29.csv'
    )))
    dose.files<-c(dose.files,paste0(path,c('dose_response_vals_all_qc.csv',
                                           'dose_response_vals_2021_05_10_all_phase_III_morpho.csv',
                                           'dose_response_vals_2021_04_26.csv',
                                           'dose_response_vals_2021_01_29.csv')))
    curv.files<-c(curv.files,paste0(path,c('fit_vals_all_qc.csv',
                                           'fit_vals_2021_05_10_all_phase_III_morpho.csv',
                                           'fit_vals_2021_04_26.csv',
                                           'fit_vals_2021_01_29.csv'
    )))
  }
  
  for(ext in extract_dirs){
    path=paste0(data.dir,'/',ext,'/')
    e.bmd<-c(e.bmd,paste0(path,'bmd_vals_2021_08_18.csv'))
    e.dose<-c(e.dose,paste0(path,'dose_response_vals_2021_08_18.csv'))
    e.curve<-c(e.curve,paste0(path,'fit_vals_2021_08_18.csv'))
  }
  
  
  ## now we read read in files
  if(length(chem.files)==3){
    #path=paste0(data.dir,'/',chem,'/')
    ##TODO: fix this to use a named list
    bmd.files<-c(bmd.files,chem.files[1])#paste0(path,'bmd_vals_all_qc.csv'))
    dose.files<-c(dose.files,chem.files[3])#paste0(path)'dose_response_vals_all_qc.csv'
    curv.files<-c(curv.files,chem.files[2])#paste0(path,'fit_vals_all_qc.csv'))
  }else{
    message("Not adding any chemical files")
  }
  #for(ext in extract.dirs){
  if(length(extract.files)==3){
    #        path=paste0(data.dir,'/',ext,'/')
    ##TODO fix this to add a named list
    e.bmd<-c(e.bmd,extract.files[1])#paste0(path,'bmd_vals_all_qc.csv'))
    e.dose<-c(e.dose,extract.files[3])##paste0(path,'dose_response_vals_all_qc.csv'))
    e.curve<-c(e.curve,extract.files[2])#paste0(path,'fit_vals_all_qc.csv'))
  }else{
    message("Not adding any environmental files")
  }
  
  chemMeta<-getChemMetadata(data.dir)
  
  sampChem<-buildSampleData(data.dir,chemMeta)
  
  endpointDetails<-getEndpointMetadata(data.dir)%>%unique()
  
  message('Processing extract response data')
  ebmds<-combineV2ChemicalEndpointData(e.bmd,is_extract=TRUE,sampChem,endpointDetails)%>%
   # select(envSampSumOutput)%>%
    unique()
  ecurves <- combineChemicalFitData(e.curve,is_extract=TRUE, sampChem,endpointDetails)%>%
    unique()
  edrs<- combineChemicalDoseData(e.dose,is_extract=TRUE, sampChem,endpointDetails)%>%
    unique()
  
  message('Processing chemical response data')
  bmds<-combineV2ChemicalEndpointData(bmd.files,is_extract=FALSE,chemMeta,endpointDetails)%>%
    unique()
  curves <-combineChemicalFitData(curv.files, is_extract=FALSE, chemMeta,endpointDetails)%>%
    unique()
  doseReps <-combineChemicalDoseData(dose.files, is_extract=FALSE, chemMeta,endpointDetails)%>%
    unique()
  
  nas<-bmds$Chemical_ID[which(is.na(bmds$AUC_Norm))]
  print(length(nas))
  
  to.remove<-setdiff(nas,sampChem$Chemical_ID)
  print(length(to.remove))
  bmds<-bmds%>%subset(!Chemical_ID%in%to.remove)
  curves<-curves%>%subset(!Chemical_ID%in%to.remove)
  doseReps<-doseReps%>%subset(!Chemical_ID%in%to.remove)
  
  ##there are mismatches, so we should figure out where those exists
  missing<-list(zebrafishNoChem=setdiff(ebmds$Sample_ID,as.character(sampChem$Sample_ID)),
                chemDataNoZebrafish=setdiff(as.character(sampChem$Sample_ID),ebmds$Sample_ID))
  
  #print(missing)
  
  ##Final output for the platform team is these 9 files
  ##chemical file - double check columns
  chems <- chemMeta%>%dplyr::select(required_chem_columns)
  write.csv(chems,file=paste0(out.dir,'chemicals.csv'),quote=TRUE,row.names=F)
  #sample file - double check columns
  samps<-sampChem%>%
    select(samp_columns)%>%
    distinct()
  
  write.csv(samps,file=paste0(out.dir,'samples.csv'),quote=T,row.names=FALSE)
  
  ##bmds 
  write.csv(bmds,file=paste0(out.dir,'zebrafishChemBMDs.csv'),quote=T,row.names = FALSE)
  write.csv(ebmds,file=paste0(out.dir,'zebrafishSampBMDs.csv'),row.names=FALSE, quote = TRUE)
  
  ##xy coords
  write.csv(curves,file=paste0(out.dir,'zebrafishChemXYCoords.csv'),row.names = FALSE, quote = TRUE)
  write.csv(ecurves,file=paste0(out.dir,'zebrafishSampXYCoords.csv'),row.names = FALSE, quote = TRUE)
  
  ##dose response
  write.csv(doseReps,file=paste0(out.dir,'zebrafishChemDoseResponse.csv'),row.names = FALSE, quote = TRUE)
  write.csv(edrs,file=paste0(out.dir,'zebrafishSampDoseResponse.csv'),row.names = FALSE, quote = TRUE)
  
  
  ##chemical to sample
  chemSamp<-sampChem%>%
    dplyr::select(sample_chem_columns)%>%
    distinct()
  write.csv(chemSamp,file=paste0(out.dir,'sampleToChemicals.csv'),row.names=FALSE, quote = TRUE)
  

}

##moving summary stats to its own function
generateSummaryStats<-function(){
  
  ##let's get all the chemicals first
  chemMeta<-read.csv(paste0(out.dir,'chemicals.csv'))
  chemClasses<-chemMeta%>%
    dplyr::select(Chemical_ID,chemical_class)%>%
    distinct()%>%
    replace_na(list(`chemical_class`='Unclassified'))
  chemClasses$chemical_class<-sapply(chemClasses$chemical_class,function(x) gsub(';deuterated','',x))
  
  ##first let's get the total chemical counts
  total.chems<-chemClasses%>%
    group_by(chemical_class)%>%
    summarize(`total chemicals`=n_distinct(Chemical_ID))
  
  ##how many samples is that chemical in? 
  sampChem<-read.csv(paste0(out.dir,'sampleToChemicals.csv'),check.names=F)%>%
    mutate(meas=as.numeric(measurement_value))%>%
    subset(!is.na(meas))%>%
    left_join(chemClasses)%>%
    group_by(chemical_class)%>%
    summarize(`Measured in samples`=n_distinct(Chemical_ID),`Number of Samples`=n_distinct(Sample_ID))
  
  #how many zebrafish experiments have that chemical?
  bmds <-read.csv(paste0(out.dir,'zebrafishChemDoseResponse.csv'),check.names=F)%>%
    group_by(Chemical_ID)%>%
    summarize(`Zebrafish endpoints`=n_distinct(`End_Point_Name`))%>%
    left_join(chemClasses)%>%
    group_by(chemical_class)%>%
    summarize(`Evaluated in Zebrafish`=n_distinct(Chemical_ID),
              `Endpoints Measured`=sum(`Zebrafish endpoints`))
  
  

  chem.eps<-total.chems%>%
    left_join(sampChem)%>%
    left_join(bmds)%>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
  
  write.table(chem.eps,paste0(out.dir,'chemCounts_v2.tsv'),row.names=F,col.names=T,sep='\t')

  
}


#' main method
#' Parsers arguments
main<-function(){
    ##now we check for additional files
    parser <- ArgumentParser()
    parser$add_argument('-s','--samples',dest='samp_files',default="",
                        help='The subsequent files are samples')
    parser$add_argument('-c','--chemicals',dest='chem_files',default='',
                        help='The subsequent files are chemicals')


    args <- parser$parse_args()
                                          #if we are adding new data, add to additional data in repo
                                        #files that we're reading in
    chem.files<-unlist(strsplit(args$chem_files,split=','))
    extract.files<-unlist(strsplit(args$samp_files,split=','))

    print(chem.files)
    print(extract.files)
    
  buildDB(chem.files,extract.files)
  generateSummaryStats()

}


test.code<-function(){
  data.dir<<-'./data/'
  out.dir<<-'./'
  buildDB()
}

main()
