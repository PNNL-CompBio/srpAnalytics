##loadGeneExpressionDAta

library(rio)
library(dplyr)
library(tidyr)

schema<-data.frame(Project=c(),
                   cas_number=c(),
                   Conc=c(),
                   Link=c(),
                   nGenes=c(),
                   Chemical_ID=c())


data.dir<-'https://raw.githubusercontent.com/PNNL-CompBio/srpAnalytics/main/data'

tab<-rio::import(paste0(data.dir,'/ZF_gex.xlsx'),skip=1)
args = commandArgs(trailingOnly=TRUE)

chem<-readr::read_csv(args[1])|>
  dplyr::select(Chemical_ID,cas_number)|>distinct()

ind_res<-tab%>% tidyr::pivot_longer(ends_with("Indication"),
                                    names_to='condition',values_to='indication')%>%
  dplyr::select(GeneID,Gene,condition,indication)%>%
  mutate(condition=stringr::str_remove(condition,'_DEG_Indication'))

##gene info - need to get basic info about gene and link to zfin db
geneinfo<-readr::read_csv(paste0(data.dir,'/allianceGenomeInfo.csv'),col_names=c('zfinId','secondId','symbol','name','organism','description'))


fc_res<-tab%>%
  tidyr::pivot_longer(ends_with("Log2FoldChange"),names_to='condition',values_to='Log2FoldChange')%>%
  dplyr::select(GeneID,Gene,condition,Log2FoldChange)%>%
  mutate(condition=stringr::str_remove(condition,'_Log2FoldChange'))

pval_res<-tab%>%
  tidyr::pivot_longer(ends_with("Adj_p_value"),names_to='condition',values_to='adj_p_value')%>%
  dplyr::select(GeneID,Gene,condition,adj_p_value)%>%
  mutate(condition=stringr::str_remove(condition,'_Adj_p_value'))

full_tab<-ind_res%>%
  full_join(fc_res)%>%
  full_join(pval_res)

##went into paper, got them from this table
##https://www.mdpi.com/1422-0067/20/10/2570

##one concentration is in ng/mL
## it is for 1746-01-6 (chem 283) - and mW is 322 g/mol
## so 1ng/mL = 1000 ug/mL 322 ~ 3 uM
## ( µg/mL ) = ( µM ) * ( MW in KD)

treat_names=data.frame(treatment=c('Acenaph_16PAH','BbF_16PAH', 'BjF_16PAH','BkF_16PAH',
                                   'DBahP_16PAH','DBaiP_16PAH','Fluoran_16PAH','Phenan_16PAH',
                                   'Retene_16PAH', 'Anthrac_16PAH','Carbazole_16PAH','X3.NF_16PAH',
                                   'X4H.CPdefP_16PAH','X9.MA_16PAH','X1.uM.Benzo.a.pyrene_PE50_AK','X10.uM.Benzo.a.pyrene_PE50_AK',
                                   'X10.uM.Dibenzo.a.l..pyrene_PE50_AK','X1.2.uM.Phenanthrene.quinone_PE50_BG',
                                   'X10.uM.Benz.a.anthracene.7.12.dione_PE50_BG','X10.uM.Benzanthrone_PE50_BG',
                                   'TCDD_TCDD','BDE.47_FRC','IPP_FRC','TBBPA.DBPE_FRC',
                                   'TBBPA_FRC','TBPH_FRC','TCEP_FRC','TCPP_FRC',
                                   'TDBPP_FRC','TPP_FRC','TiBP_FRC'),
                       Conc=c('50uM','50uM','14.9uM','1.9uM','10uM','5uM','50uM','12.2uM',
                              '50uM','50uM','1.9uM','16.2uM','50uM','50uM','1uM','10uM',
                              '10uM','2uM','10uM','10uM','3uM','85uM','19.8uM','4uM',
                              '4uM','72uM','85uM','85uM','3uM','8uM','158uM'),
                       cas_number=c('83-32-9','205-99-2','205-82-3','207-8-9','189-64-0','189-55-9','206-44-0','65996-93-2',
                                    '483-65-8','120-12-7','86-74-8','892-21-7','5737-13-3','779-02-2','50-32-8','50-32-8',
                                    '191-30-0','84-11-7','2498-66-0','82-05-3','1746-01-6','5436-43-1','68937-41-7','21850-44-2',
                                    '79-94-7','26040-51-7','115-96-8','13674-84-5','126-72-7','115-86-6','126-71-6'))


##differential expression counts for condition

allgenes<-full_tab|>
  tidyr::separate(condition,into=c('control','treatment'),sep='_v_')|>
  left_join(treat_names)|>
  dplyr::select(-c(control,treatment))|>
  mutate(Project='Zebrafish',Link='')|>
  left_join(chem)|>
  mutate(Concentration=as.numeric(stringr::str_replace(Conc,'uM','')))|>
  select(-Conc)


diffex <- allgenes|>
  subset(indication==1)

res<-allgenes|>
  subset(indication!=0)|>
  mutate(up=Log2FoldChange>0)|>
  group_by(Project,cas_number,Concentration,Link,Chemical_ID,up)|>
  summarize(nGenes=n())|>subset(!is.na(up))


res<-res|>
  tidyr::pivot_wider(names_from='up',
                     values_from='nGenes')|>
  rename(UpRegulatedGenes='TRUE',DownRegulatedGenes='FALSE')


##need to get mapping to drug name
write.table(res,file ='/tmp/srpDEGstats.csv',sep=',',quote=F,row.names=F,col.names=T)
write.table(allgenes,file='/tmp/allGeneEx.csv',sep=',',quote=F,row.names=F,col.names=T)


