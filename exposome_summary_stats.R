##get exposome data

library(jsonlite)
library(httr)



portal_name='https://montilab.bu.edu/Xposome-API/portals'



#test with a chemical id - I'm assuming it's CAS ids
chemical_id='524-42-5'

#' get GO terms for each chemical id
getGoTerms<-function(chemical_id){

  ##get gsva stats
  url2 <- paste0("https://montilab.bu.edu/Xposome-API/gs_enrichment?portal=", 
               portal_name, "&chemical_id=", chemical_id, 
               "&geneset=Hallmark&gsva=gsva&summarize.func=median")

  # Send GET Request to API
  res <- GET(url = url2, encode = 'json')

  stop_for_status(res)

  gs_enrichment_stat <- fromJSON(fromJSON(rawToChar(res$content)))

}


#getdiffex genes
getGenes<-function(chemical_id){

  url1 <- paste0("https://montilab.bu.edu/Xposome-API/gene_expression?portal=", 
                 portal_name, "&chemical_id=", chemical_id, 
                 "&summarize.func=median&landmark=TRUE&do.markers=TRUE&do.scorecutoff=TRUE")
  
  # Send GET Request to API
  res <- GET(url = url1, encode = 'json')
  
  stop_for_status(res)
  
  gene_expression_stat <- fromJSON(fromJSON(rawToChar(res$content)))
  
  
}