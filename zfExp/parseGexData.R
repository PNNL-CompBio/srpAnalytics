##loadGeneExpressionDAta

library(rio)
library(dplyr)
library(tidyr)
library(enrichR)
library(readr)

schema<-data.frame(Project=c(),
                   cas_number=c(),
                   Conc=c(),
                   Link=c(),
                   nGenes=c(),
                   Chemical_ID=c())


data.dir<-'https://raw.githubusercontent.com/PNNL-CompBio/srpAnalytics/main/data/zfExp'
out.dir='/tmp/'

generateGeneExamples<-function(genelist,chems){
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  genelist<-genelist|>
    subset(!is.na(Gene))|>
    mutate(Significant=adj_p_value<0.05,LogPValue=-log10(adj_p_value),Direction=ifelse(Log2FoldChange>0,'Up-regulated','Down-regulated'))

  ccounts<-genelist|>left_join(chems)|>
    dplyr::select(cas_number,chemical_class)|>
    distinct()|>
    group_by(chemical_class)|>
    summarize(chem_counts=n())|>
    mutate(chem_frac=chem_counts/sum(chem_counts))
  ccounts$ymax=cumsum(ccounts$chem_frac)
  ccounts$ymin=c(0, head(ccounts$ymax, n=-1))

  ggplot(ccounts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=chemical_class)) +
      geom_rect() +
      coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
      xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
    theme_void()
  ##need gene summar plots
  testgenes<-sample(genelist$Gene,5)

  for(gene in testgenes){
    p1<-subset(genelist,Gene==gene)|>ggplot(aes(x=Log2FoldChange,y=LogPValue,col=Significant,label=cas_number))+geom_label_repel()+geom_point()
    p2<-subset(genelist,Gene==gene)|>ggplot(aes(x=Direction,fill=Significant))+geom_bar(position='dodge')
   # tab<-subset(genelist,Gene==gene)|>dplyr::select(cas_number,Chemical_ID,Concentration,Log2FoldChange,adj_p_value)|>grid.table()
    cowplot::plot_grid(cowplot::plot_grid(p1,p2,nrow=1))#,tab,nrow=2)
    ggsave(paste0('genePlotsFor_',gene,'.png'),width=14)

  }


}

#test script to try out some examples
generateChemicalExamples<-function(genelist,deglist){
  library(ggplot)
  library(cowplot)
  testchems<-sample(pathlist$Chemical_ID,5)

  for(chem in testchems){
    genes<-subset(genelist,Chemical_ID==chem)|>
      mutate(Pvalue=-1*log10(adj_p_value),Significant=ifelse(indication==1,TRUE,FALSE),Concentration=as.factor(Concentration))|>
      ggplot(aes(x=Log2FoldChange,y=Pvalue,col=Significant,alpha=0.5,shape=Concentration))+geom_point()

    genecount<-subset(deglist,Chemical_ID==chem)|>
      mutate(Direction=ifelse(Log2FoldChange>0,"Up","Down"),Significant=ifelse(indication==1,TRUE,FALSE),Concentration=as.factor(Concentration))|>
      subset(Significant)|>
      subset(!is.na(Direction))|>
      ggplot(aes(x=Direction,fill=Concentration))+geom_bar(position='dodge')

    paths<-subset(deglist,Chemical_ID==chem)|>
      subset(toPlot==1)|>
      tidyr::separate(Overlap,into=c('top','bottom'),sep='/')|>
      mutate(geneCount=as.numeric(top)/as.numeric(bottom),Concentration=as.factor(Concentration))|>
      ggplot(aes(x=reorder(Term,geneCount),y=geneCount,fill=-1*log10(Adjusted.P.value)))+geom_bar(stat='identity',position='dodge')+coord_flip()
    ##volcano plot for genes
    ##barplot for pathways

    cowplot::plot_grid(cowplot::plot_grid(genes,genecount,nrow=1),paths,nrow=2)
    ggsave(paste0('chemPlotsFor_',chem,'.png'))

  }

}

##we only can plot the top20 or so, so including a value here to enable filtering on web end
enrichSelectTop<-function(genelist,path,pvalue=0.05,top=20){

  res<-enrichr(genelist,path)
  res<-res[[path]]

  res<-res|>
    subset(Adjusted.P.value<pvalue)|>
    arrange(Adjusted.P.value)

  res$toPlot<-rep(1,nrow(res))
  if(nrow(res)>top){
    res$toPlot<-rep(0,nrow(res))
    res$toPlot[1:top]<-1
  }
  return(res)
}
##call functional enrichment and store results
doEnrich<-function(genelist){
    condition<-genelist|>dplyr::select(Chemical_ID,Concentration)|>
      distinct()
    setEnrichrSite("FishEnrichr")
#    dbs<-listEnrichrDbs()$libraryName
    path=c('Coexpression_Predicted_GO_Biological_Process_2018')

    allpaths<-genelist|>
      subset(indication==1)|>
    #  subset(Chemical_ID%in%c(3138,3130,3148))|>
      group_by(Chemical_ID,Concentration)|>
      summarize(enrich=enrichSelectTop(Gene,path,0.05,20))

    ##now unnest and filter
    sigpaths<-allpaths|>
      unnest(cols=c(enrich))

    ##filter for signifiance, then move to long form table
    return(sigpaths)
}


##main function to do all the things
main<-function(args=c()){
  if(length(args)!=3)
    args = commandArgs(trailingOnly=TRUE)
  #print(args)
  if(length(args)<3){
      print('Need to call script with path to GEX files (comma delimited) and chemicals.csv and gene info file')
      quit()
  }


  tab<-rio::import(args[1],which=1,skip=1)###so far only equipped to handle one gene expression file

  chem<-readr::read_csv(args[2])|>
    dplyr::select(Chemical_ID,cas_number)|>distinct()

  ind_res<-tab%>% tidyr::pivot_longer(ends_with("Indication"),
                                      names_to='condition',values_to='indication')%>%
    dplyr::select(GeneID,Gene,condition,indication)%>%
    mutate(condition=stringr::str_remove(condition,'_DEG_Indication'))

  ##gene info - need to get basic info about gene and link to zfin db
  geneinfo<-readr::read_csv(args[3],col_names=c('zfinId','secondId','symbol','name','organism','description'))


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


  enrich<-doEnrich(allgenes)

  res<-res|>left_join(enrich)
  ##need to get mapping to drug name
  readr::write_csv(enrich,file =paste0(out.dir,'srpDEGPathways.csv'))
  readr::write_csv(res,file =paste0(out.dir,'srpDEGstats.csv'))
  readr::write_csv(allgenes,file=paste0(out.dir,'allGeneEx.csv'))
  return(list(deg=res,genes=allgenes))
}



main()#args<-c('../data/zfExp/ZF_gex.xlsx','../data/chemicalIdMapping.csv','../data/zfExp/allianceGenomeInfo.csv'))

