##get exposome data

library(jsonlite)
library(httr)
library(dplyr)
library(tidyr)
#library(cowplot)
url <- "https://montilab.bu.edu/Xposome-API/projects?all=Yes"
res <- GET(url = url, encode = 'json')
stop_for_status(res)

projects <- fromJSON(fromJSON(rawToChar(res$content)))

print(paste('We now have data from',length(projects),'projects'))

#portal_name=  'https://montilab.bu.edu/Xposome-API/portals'

##read in all chemicals
all.chems <- read.table('chemicalIdMapping.csv',sep=',',header=T,fileEncoding = "UTF-8-BOM")


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
      gs_enrichment_stat <- data.frame(`GenesetName`=NULL,`Geneset`=NULL,
                                       `Summary Score`=NULL,`GSScore`=NULL,`Conc`=NULL)
  else{
    gs_enrichment_stat <- fromJSON(fromJSON(rawToChar(res$content)))
    gs_enrichment_stat <- gs_enrichment_stat%>%
      tibble::rownames_to_column('GenesetName')
    gs_enrichment_stat <- gs_enrichment_stat%>%
      pivot_longer(cols=grep('GS Score',colnames(gs_enrichment_stat)),names_to='Conc',values_to='GSScore')%>%
      mutate(Conc=stringr::str_replace_all(Conc,'GS Score ',''))

    }
  print(head(gs_enrichment_stat))
  data.frame(Project=proj,cas_number=chemical_id,gs_enrichment_stat)
}


#' getGenes
#' Gets list of genes differentially expressed for a particular chemical and project
#' @param chemical_id CAS id
#' @param proj project id
#' @return data frame
getGenes<-function(chemical_id,proj){

  url1 <- paste0("https://montilab.bu.edu/Xposome-API/gene_expression?project=",
                 proj, "&chemical_id=", chemical_id,"&landmark=FALSE&do.scorecutoff=FALSE")
  # Send GET Request to API
  res <- GET(url = url1, encode = 'json')
  chem_gene_link=paste0("https://montilab.bu.edu/Xposome/?page=",proj,
                   "&tab=chemical_explorer&chemical_id=",chemical_id,"&stat=gene_expression")

  if(res$status_code!=200)
    gene_expression_stat=data.frame(Gene=NULL,GeneName=NULL,Direction=NULL,`Summary Score`=NULL,`Zscore`=NULL,`Conc`=NULL,Link=NULL)
  else{
    gene_expression_stat <- fromJSON(fromJSON(rawToChar(res$content)))
    gene_expression_stat <- gene_expression_stat%>%
      tibble::rownames_to_column('GeneName')
    gene_expression_stat <-gene_expression_stat%>%
      pivot_longer(cols=grep('ModZScore',colnames(gene_expression_stat)),names_to='Conc',values_to='ModZScore')%>%
      mutate(Conc=stringr::str_replace_all(Conc,'ModZScore ',''))
    if('Landmark_Gene'%in%colnames(gene_expression_stat))
      gene_expression_stat <- select(gene_expression_stat,-Landmark_Gene)
  }
 # print(head(gene_expression_stat))
    data.frame(Project=proj,cas_number=chemical_id,gene_expression_stat,Link=chem_gene_link)

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
  print(paste('Found',length(overlap),'cas ids in common in project',proj))

   #chem_Term_link=paste0("https://montilab.bu.edu/Xposome/?page=",proj,
  #                 "&tab=chemical_explorer&chemical_id=",chem,"&stat=gene_set_enrichment")

  gg=do.call(rbind,lapply(overlap,function(chem) getGenes(chem,proj)))
#  gt=do.call(rbind,lapply(overlap,function(chem) getGoTerms(chem,proj)))
  full.list[[proj]]<-list(genes=gg)#,goterms=gt)
}

gene.tab<-NULL
term.tab<-NULL
for(i in 1:length(projects)){
  print(i)
  gene.tab <- rbind(gene.tab, full.list[[i]]$genes)
 # term.tab <-rbind(term.tab,full.list[[i]]$goterms)
}


sig.genes <- gene.tab%>%rowwise()%>%mutate(absVal=abs(ModZScore))%>%
  subset(absVal>1.63)

## library(ggplot2)
## chems<-unique(gene.tab$cas_number)[100:105]
## for(i in chems){
##   p1<-sig.genes%>%
##     subset(cas_number==i)%>%
##     ggplot()+geom_jitter(aes(x=Conc,y=ModZScore,col=Project,alpha=0.5))+ggtitle(i)+ scale_x_discrete(guide = guide_axis(angle = 90))
##   p2<-sig.genes%>%
##     subset(cas_number==i)%>%
##     ggplot()+geom_bar(aes(x=Conc,fill=Project),position='dodge')+ggtitle(i) +scale_x_discrete(guide = guide_axis(angle = 90))

##   res=cowplot::plot_grid(plotlist=list(p1,p2))
##   ggsave(paste0('summary',i,'.png'),plot=res)
## }


map <-all.chems%>%
  dplyr::select(cas_number,Chemical_ID)%>%distinct()


sg.stats <- sig.genes%>%
  group_by(Project,cas_number,Conc,Link)%>%
  summarize(nGenes=n_distinct(Gene))%>%
  left_join(map)

write.table(sg.stats,file='/tmp/sigGeneStats.csv',sep=',',row.names=F)
##not using this for now:
#write.table(sig.genes,file='data/sigGeneExp.csv',sep=',',row.names=F)
