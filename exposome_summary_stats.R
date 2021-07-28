##get exposome data

library(jsonlite)
library(httr)

url <- "https://montilab.bu.edu/Xposome-API/projects?all=Yes"
res <- GET(url = url, encode = 'json')
stop_for_status(res)

projects <- fromJSON(fromJSON(rawToChar(res$content)))

print(paste('We now have',length(projects),'projects'))

#portal_name=  'https://montilab.bu.edu/Xposome-API/portals'

##read in all chemicals
all.chems <- read.table('data/chemicalIdMapping.csv',sep=',',header=T,fileEncoding = "UTF-8-BOM")


#' get GO terms for each chemical id
getGoTerms<-function(chemical_id,proj){

  ##get gsva stats
  url2 <- paste0("https://montilab.bu.edu/Xposome-API/gs_enrichment?project=", 
               proj, "&chemical_id=", chemical_id)#, 
               #"&geneset=Hallmark&gsva=gsva&summarize.func=median")

  print(url2)
  # Send GET Request to API
  res <- GET(url = url2, encode = 'json')

  
  if(res$status_code!=200)
    return(NULL)
  
  gs_enrichment_stat <- fromJSON(fromJSON(rawToChar(res$content)))
  gs_enrichment_stat
}


#getdiffex genes
getGenes<-function(chemical_id,proj){

  url1 <- paste0("https://montilab.bu.edu/Xposome-API/gene_expression?project=", 
                 proj, "&chemical_id=", chemical_id)#, 
                 #"&summarize.func=median&landmark=TRUE&do.markers=TRUE&do.scorecutoff=TRUE")
  
  print(url1)
  # Send GET Request to API
  res <- GET(url = url1, encode = 'json')
  
  if(res$status_code!=200)
    return(NULL)res
  
  gene_expression_stat <- fromJSON(fromJSON(rawToChar(res$content)))
  gene_expression_stat 
  
}

full.list <- list()
for(proj in projects){
  # Now we get the list of chemicals available
  chem_list_url <- paste0("https://montilab.bu.edu/Xposome-API/chemicals?projects=",proj,"&chemical_ids=all")
  # Send GET Request to API
  res <- GET(url = chem_list_url, encode = 'json')
  if(res$status_code=='200')
    chemicals <- fromJSON(fromJSON(rawToChar(res$content)))
  
  overlap <- intersect(all.chems$cas_number,chemicals$CAS)
  print(paste('Found',length(overlap),'cas ids in common'))
  
  
  gg=getGenes(paste(overlap,collapse=','),proj)
  gt=getGoTerms(paste(overlap,collapse=','),proj)
  full.list[[proj]]<-list(genes=gg,goterms=gt)
}