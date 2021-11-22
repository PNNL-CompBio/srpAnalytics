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
#data.dir='./data/'
required_sample_columns<-c("SampleNumber","date_sampled","sample_matrix","technology",
                           "Sample_ID","ProjectName","SampleName","LocationLat",
                           "LocationLon","LocationName","LocationAlternateDescription",
                           "AlternateName")

required_chem_columns<-c()

required_mapping_columns<-c("SampleNumber","Sample_ID",'zf_lims_id')## added these to complete mapping
#extra columns: "bioassay_sample_from_lims","parent_samples_from_lims","child_samples_from_lims")

required_bmd_columns<-list(bmd=c('Chemical_ID','End_Point','Model','BMD10','BMD50',"Min_Dose","Max_Dose",
                                "AUC_Norm","DataQC_Flag","BMD_Analysis_Flag"),#,"BMD10_Flag","BMD50_Flag"),
                          doseRep=c("Chemical_ID","End_Point","Dose","Response","CI_Lo","CI_Hi"),
                          fitVals=c("Chemical_ID","End_Point","X_vals","Y_vals"))

##output directory is fixed
out.dir<-'/tmp/'
#out.dir<-'./'
#' Get chemical metadata, which is stored in `data.dir`
#' @param data.dir path to standardized data that enables matching across datasets
#' @return data.frame
getChemMetadata<-function(data.dir){
    ##This mapping file consumes data from OSU to match identifiers to CAS
    fullMapping<-read.csv2(paste0(data.dir,'chemicalIdMapping.csv'),header=TRUE,sep=',',fileEncoding="UTF-8-BOM")%>%
        dplyr::select(cas_number,zf.cid,Chemical_ID,chemical_class)%>%
        distinct()

    curatedDesc <- readxl::read_xlsx(paste0(data.dir,'ChemicalDescriptions.xlsx'))%>%
        select(cas_number='CASRN',chemDescription='USE CATEGORY/DESCRIPTION')

    ##this file comes from COMPTOX
    ##here we join the chemical metadata from the comptox dashboard
    chemMeta<-readxl::read_xls(paste0(data.dir,
                                      'CompToxChemicalsDashboard-Batch-Search_2021-11-05_17_57_46.xls'))%>%
        select(INPUT,DTXSID,PREFERRED_NAME,INCHIKEY,SMILES,MOLECULAR_FORMULA,
               AVERAGE_MASS,PUBCHEM_DATA_SOURCES)%>%
        rename(cas_number='INPUT')%>%
        full_join(fullMapping)%>% ## add in ID mapping information -- THIS WILL INTRODUCE NAs
        left_join(getNewChemicalClass(data.dir))%>% ##add in chemical class
        left_join(curatedDesc)%>% ## add in description data
        tidyr::replace_na(list(newClass='Unclassified'))%>% ##allow for unclassified
        tidyr::replace_na(list(PREFERRED_NAME='Chemical name unknown'))%>% ##should we call this something else?
        mutate(PREFERRED_NAME=stringr::str_replace_all(PREFERRED_NAME,'^-$',"Chemical name unknown"))%>%
        select(-c(chemical_class,Classification))%>%
        rename(chemical_class='newClass')

    nocas=grep("NOCAS",chemMeta$cas_number)
    if(length(nocas)>0){
        message(paste0('removing ',length(nocas),'chems with no cas'))
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
        subset(`Portal Display`=='Display')%>%
        rename(End_Point='Abbreviation',`End Point Name`='Simple name (<20char)')%>%
        rename(endPointLink='Ontology Link')%>%
        select(-`Portal Display`)%>%
      mutate(End_Point=stringr::str_trim(End_Point))

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
    sampChem <- read.csv(paste0(data.dir,'/fses_data_for_pnnl_4-27-2021.csv'),fileEncoding="UTF-8-BOM")%>%
                                        #  sampChem<-read.csv(paste0(data.dir,'/pnnl_bioassay_sample_query_1-14-2021.csv'))%>%
        subset(SampleNumber!='None')%>%
        subset(cas_number!='NULL')%>%
        mutate(water_concentration_molar=stringr::str_replace_all(water_concentration_molar,'BLOD|NULL|nc:BDL',"0"))%>%
        mutate(measurement_value_molar=stringr::str_replace_all(measurement_value_molar,'BLOD|NULL|BDL',"0"))%>%
        mutate(water_concentration=stringr::str_replace_all(water_concentration,'BLOD|NULL|BDL',"0"))%>%
                                        # subset(water_concentration_molar!='0.0')%>%
        subset(!measurement_value_molar%in%c('0'))%>%
        subset(!measurement_value%in%c("0","NULL",""))%>%
        select(-c(Sample_ID,Chemical_ID)) ##These two are added in the 4/27 version of the file

    chemDat<-chemMeta%>%
        select(Chemical_ID,cas_number,AVERAGE_MASS)%>%
        distinct()

    ##This mapipng file maps tanguay lab identifiers to those in the anderson lab
    ids<-read.csv(paste0(data.dir,'/sampleIdMapping.csv'),fileEncoding="UTF-8-BOM")%>%
        select(Sample_ID,SampleNumber)%>%
        distinct()

    sampChem <-sampChem%>%
        inner_join(chemDat,by='cas_number')%>%
        left_join(ids,by='SampleNumber')%>%
        distinct()

    ##now we have to fix duplicate sample names

    #get duplicates, assign
    all.samp.names<-sampChem%>%
      select(Sample_ID,SampleName)%>%
      distinct()%>%
      mutate(isDupe=duplicated(SampleName))

    #filter only duplicates
    dupe.samp.names<-all.samp.names%>%
      subset(isDupe)%>%
      arrange(SampleName)
    #get number
    num.dupes<-dupe.samp.names%>%group_by(SampleName)%>%
      summarize(nid=n_distinct(Sample_ID))%>%left_join(dupe.samp.names)
    #paste sequence to end of each name
    new.names<-num.dupes%>%
      select(SampleName,nid)%>%distinct()%>%
      rowwise()%>%
      mutate(newName=paste(SampleName,seq(nid),collapse=','))%>%
      tidyr::separate_rows(newName,sep=',')%>%
      select(oldSN='SampleName',newName)%>%
      arrange(oldSN)
    #combine in new data frame
    full.rep=data.frame(dupe.samp.names,new.names)%>%
      select(Sample_ID,SampleName='newName')

    new.samp.names<-all.samp.names%>%subset(!isDupe)%>%
      select(-isDupe)%>%
      rbind(full.rep)

    sampChem<-sampChem%>%
      select(-SampleName)%>%#remove duplicates
      left_join(new.samp.names)

    print("Updating concentration data")
                                        #TODO: separate those without molar concentrations, recompute
    blanks<-sampChem%>%subset(measurement_value_molar=='')
    nonblanks<-sampChem%>%subset(measurement_value_molar!='')
    fixed<-blanks%>%
        mutate(concentration=as.numeric(measurement_value))%>%
        mutate(AVERAGE_MASS=as.numeric(AVERAGE_MASS))%>%
        rowwise()%>%
        mutate(measurement_value_molar=concentration*1000/AVERAGE_MASS)%>%
        select(-concentration)

    sampChem<-rbind(nonblanks,fixed)%>%
        select(-AVERAGE_MASS)%>%
        distinct()

    ##now we have one more rename of samples and metadata
    sampleNameRemap<-readxl::read_xlsx(paste0(data.dir,'/envSampCleanMapping.xlsx'))%>%
      dplyr::select(Sample_ID,date_sampled,sample_matrix,technology,ProjectName,SampleName='NewSampleName',
                    LocationName='NewLocationName')%>%
      distinct()

    sampChem<-sampChem%>%
      select(-c(sample_matrix,technology,SampleName,LocationName,date_sampled))%>% ### remove old names
      left_join(sampleNameRemap)%>% ##add in updates
      distinct()

    return(sampChem)

}
#'combineChemicalEndpointData produces the summary statistics from the BMD analysis
#'@param bmdfiles a list of files that come from the BMD pipeline
#'@param is_extract - if our chemical id is an extract we use a different metadata
#'@return a data.frame
##We will release an 'endpoint file for each condition'
combineChemicalEndpointData<-function(bmdfiles,is_extract=FALSE,sampChem,endpointDetails){

  print(paste('Combining bmd files:',paste(bmdfiles,collapse=',')))
  cols <- required_bmd_columns$bmd
  files <- lapply(bmdfiles,function(x) read.csv(x)%>%dplyr::select(cols))

  mid.bmd<-do.call(rbind,files)%>%
      dplyr::select(cols)


  dupes<-which(mid.bmd%>%select(Chemical_ID,End_Point)%>%duplicated())
  if(length(dupes)>0){
    mid.bmd<-mid.bmd[-dupes,]
  }
  if(is_extract){
    sdSamp<-sampChem%>%tidyr::separate('Sample_ID',into=c('tmpId','sub'),sep='-',remove=FALSE)%>%
      select(-sub)

    full.bmd<-mid.bmd%>%
      dplyr::mutate(tmpId=as.character(Chemical_ID))%>%
      dplyr::select(-Chemical_ID)%>%
      full_join(sdSamp,by='tmpId')#%>%#%>%mutate(Chemical_ID<-as.character(zaap_cid)))%>%

    #fix up sample ids
    nas<-which(is.na(full.bmd$Sample_ID))
    full.bmd$Sample_ID[nas]<-full.bmd$tmpId[nas]

    #now fix up sample names
    new.nas<-which(is.na(full.bmd$SampleName))
    full.bmd$SampleName[new.nas]<-paste('Sample',full.bmd$Sample_ID[new.nas])

    full.bmd<-full.bmd%>%
      right_join(endpointDetails)%>%
      dplyr::select(-c('End_Point','tmpId'))%>%
    distinct()%>%
      tidyr::replace_na(list(LocationName='None'))
}
  else{
    full.bmd<-mid.bmd%>%
      #dplyr::mutate(`Chemical_ID`=as.character(Chemical_ID))%>%
      inner_join(sampChem)%>%
#      rename(Chemical_ID<-'zf.cid')%>%
      right_join(endpointDetails)%>%
      distinct()%>%select(-c('End_Point'))%>%
      tidyr::replace_na(list(chemical_class='Unclassified'))##should we remove endpoint YES
  }

  ##now we fix QC values
  full.bmd <- full.bmd%>%
    rename(qc_num='DataQC_Flag')%>%
    mutate(DataQC_Flag=ifelse(qc_num%in%c(0,1),'Poor',ifelse(qc_num==4,'Moderate','Good')))%>%
    rowwise()%>%
    mutate(Model=stringr::str_replace_all(Model,"NULL","None"))%>%
      select(-c(qc_num,BMD_Analysis_Flag))

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
        mutate(zf.cid=as.character(Chemical_ID))%>%
        rename(ChemicalId='zf.cid')%>%
        subset(X_vals!="NULL")%>%
        mutate(X_vals=as.numeric(X_vals))%>%
        mutate(Y_vals=as.numeric(Y_vals))

    })

    mid.bmd<-do.call(rbind,fixed.files)

   # mid.bmd<-do.call(rbind,files)%>%
    #    dplyr::select(required_bmd_columns$fitVals)%>%
    #    mutate(zf.cid=as.character(Chemical_ID))%>%
    #    rename(ChemicalId='zf.cid')%>%
    #    subset(X_vals!="NULL")%>%
    #    mutate(X_vals=as.numeric(X_vals))%>%
    #    mutate(Y_vals=as.numeric(Y_vals))

    ##this alone fails to determine which are the duplicates...
    #dupes<-which(mid.bmd%>%
    #             select(Chemical_ID,End_Point,X_vals)%>%duplicated())

    #if(length(dupes)>0){
       # print(mid.bmd[dupes,])
    #  mid.bmd<-mid.bmd[-dupes,]
    #}
    print(head(mid.bmd))
    full.bmd<-mid.bmd%>%
        right_join(endpointDetails)%>%
        distinct()%>%
        select(-c(End_Point,Description,ChemicalId))
    print(head(full.bmd))
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
         mutate(zf.cid=as.character(Chemical_ID))%>%
       rename(ChemicalId='zf.cid')
    })

    mid.bmd<-do.call(rbind,fixed.files)
    # mid.bmd<-do.call(rbind,files)%>%
    #     dplyr::select(required_bmd_columns$doseRep)%>%
    #     mutate(zf.cid=as.character(Chemical_ID))%>%
    #     rename(ChemicalId='zf.cid')
    #
    #
    # dupes<-which(mid.bmd%>%select(Chemical_ID,End_Point,Dose)%>%
    #              mutate(Dose=as.numeric(Dose))%>%
    #              duplicated())
    #
    # if(length(dupes)>0){
    #   mid.bmd<-mid.bmd[-dupes,]
    # }
    #
    #
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

    }else{
        full.bmd<-full.bmd%>%
        subset(Chemical_ID%in%sampChem$Chemical_ID)
       # full.bmd <- rename(full.bmd,Sample_ID='Chemical_ID')
   }
  return(unique(full.bmd))
}


chem_dirs=c('phase_I_II_III')#,'zf_morphology')
extract_dirs=c('extracts')


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

    bmd.files<-c()
    dose.files<-c()
    curv.files<-c()

    e.bmd<-c()
    e.dose<-c()
    e.curve<-c()

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
        print("Not adding any environmental files")
    }

    chemMeta<-getChemMetadata(data.dir)
    sampChem<-buildSampleData(data.dir,chemMeta)

    sampData<-sampChem%>%
        dplyr::select(required_sample_columns)%>%
        distinct()

    sampToChem<-sampChem%>%
        dplyr::select(Sample_ID,Chemical_ID,SampleNumber,SampleName,cas_number,
                      measurement_value_molar,measurement_value_molar_unit,
                      water_concentration,water_concentration_unit,
                      water_concentration_molar,water_concentration_molar_unit)%>%
        distinct()


    print('Processing extract response data')
    endpointDetails<-getEndpointMetadata(data.dir)%>%unique()
    ebmds<-combineChemicalEndpointData(e.bmd,is_extract=TRUE,sampData,endpointDetails)%>%
      unique()
    ecurves <- combineChemicalFitData(e.curve,is_extract=TRUE, sampData,endpointDetails)%>%
      unique()
    edrs<- combineChemicalDoseData(e.dose,is_extract=TRUE, sampData,endpointDetails)%>%
      unique()

    print('Processing chemical response data')
    bmds<-combineChemicalEndpointData(bmd.files,is_extract=FALSE,chemMeta,endpointDetails)%>%
      unique()
    curves <-combineChemicalFitData(curv.files, is_extract=FALSE, chemMeta,endpointDetails)%>%
      unique()
    doseReps <-combineChemicalDoseData(dose.files, is_extract=FALSE, chemMeta,endpointDetails)%>%
      unique()

    nas<-bmds$Chemical_ID[which(is.na(bmds$AUC_Norm))]
    print(length(nas))
    to.remove<-setdiff(nas,sampToChem$Chemical_ID)
    print(length(to.remove))
    bmds<-bmds%>%subset(!Chemical_ID%in%to.remove)
    curves<-curves%>%subset(!Chemical_ID%in%to.remove)
    doseReps<-doseReps%>%subset(!Chemical_ID%in%to.remove)


    #print(head(curves))
### we used to delete duplicates here, but have pushed
    #that back into other functions
    #bdupes=curves%>%select('End Point Name','X_vals','Chemical_ID')%>%duplicated()
    #print(bdupes)
    #if(length(bdupes)>0)
    #    curves<-curves[-bdupes,]

    #edupes<-ecurves%>%select('End Point Name','X_vals','Sample_ID')%>%duplicated()
    #if(length(edupes)>0)
    #    ecurves<-ecurves[-edupes,]
    #print(head(curves))
    ##there are mismatches, so we should figure out where those exists
    missing<-list(zebrafishNoChem=setdiff(ebmds$Sample_ID,as.character(sampChem$Sample_ID)),
                  chemDataNoZebrafish=setdiff(as.character(sampChem$Sample_ID),ebmds$Sample_ID))

    #print(missing)

    ##Final output for the platform team is these 4 files
    write.csv(bmds,file=paste0(out.dir,'chemSummaryStats.csv'),quote=T,row.names = FALSE)
    write.csv(ebmds,file=paste0(out.dir,'envSampSummaryStats.csv'),row.names=FALSE, quote = TRUE)

    write.csv(curves,file=paste0(out.dir,'chemXYcoords.csv'),row.names = FALSE, quote = TRUE)
    write.csv(ecurves,file=paste0(out.dir,'envSampXYcoords.csv'),row.names = FALSE, quote = TRUE)

    write.csv(doseReps,file=paste0(out.dir,'chemdoseResponseVals.csv'),row.names = FALSE, quote = TRUE)
    write.csv(edrs,file=paste0(out.dir,'envSampdoseResponseVals.csv'),row.names = FALSE, quote = TRUE)

    write.csv(sampToChem,file=paste0(out.dir,'chemicalsByExtractSample.csv'),row.names=FALSE, quote = TRUE)

    ##let's do one last summary
    samp.count<-sampToChem%>%select(Sample_ID,Chemical_ID)%>%
      group_by(Chemical_ID)%>%
      summarize(`Number of samples`=n_distinct(Sample_ID))

    chem.counts<-bmds%>%
      select(c('Chemical_ID','chemical_class','End Point Name','AUC_Norm'))%>%
        subset(!is.na(AUC_Norm))%>%
        group_by(chemical_class)%>%
        summarize(`Chemicals`=n_distinct(Chemical_ID))

    chem.eps<-bmds%>%
        group_by(Chemical_ID,chemical_class)%>%
        summarize(`End Points`=n_distinct(`End Point Name`))%>%
        full_join(samp.count)%>%
        tidyr::replace_na(list(`End Points`=0,`Number of samples`=0,`chemical_class`='None'))%>%
        group_by(chemical_class)%>%
        summarize(`End Points`=sum(`End Points`),`Samples`=sum(`Number of samples`))%>%
        left_join(chem.counts)

    chem.count<-sampToChem%>%select(Sample_ID,Chemical_ID)%>%
      group_by(Sample_ID)%>%
      summarize(`Number of chemicals`=n_distinct(Chemical_ID))

    samp.counts<-ebmds%>%
      group_by(LocationName)%>%
      summarize(`Number of samples`=n_distinct(Sample_ID))

    samp.eps<-ebmds%>%
     select(c('Sample_ID','LocationName','End Point Name','AUC_Norm'))%>%
      subset(!is.na(AUC_Norm))%>%
    group_by(Sample_ID,LocationName)%>%
    summarize(`End Points`=n_distinct(`End Point Name`))%>%
    full_join(chem.count)%>%
    tidyr::replace_na(list(`End Points`=0,`Number of chemicals`=0,`LocationName`='None'))%>%
      group_by(LocationName)%>%
      summarize(`End Points`=sum(`End Points`),Chemicals=sum(`Number of chemicals`))%>%
      left_join(samp.counts)

    write.table(chem.eps,paste0(out.dir,'chemCounts.Rmd'),row.names=F,col.names=T,sep='|')
    write.table(samp.eps,paste0(out.dir,'sampCounts.Rmd'),row.names=F,col.names=T,sep='|')

}

main()
