fixDateColoumn<-function(fname){
  tab<-readr::read_csv(fname)
  datecol<-tab$date_sampled
  
  mornings<-grep("AM",datecol)
  evenings<-grep("PM",datecol)
  other<-setdiff(1:length(datecol),union(mornings,evenings))
  
  monthlist<-c('01','02','03','04','05','06','07','08','09','10','11','12')
  names(monthlist)<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov',"Dec")
  morntab<-data.frame(od=datecol[mornings])|>
    tidyr::separate(od,into=c('month','day','year','hour','min'),remove=FALSE)|>
    dplyr::mutate(day=ifelse(as.numeric(day)<10,paste0('0',day),day))|>
    dplyr::mutate(month=monthlist[month])|>
    tidyr::unite(year,month,day,col='date',sep='-')|>
    dplyr::mutate(min=gsub('AM',':00',min))|>
    dplyr::mutate(hour=ifelse(hour=='12','00',hour))|>
    tidyr::unite(hour,min,sep=':',col='time')|>
    tidyr::unite(date,time,col=newdate,sep=' ')
  
  evetab<-data.frame(od=datecol[evenings])|>
    tidyr::separate(od,into=c('month','day','year','hour','min'),remove=FALSE)|>
    dplyr::mutate(day=ifelse(as.numeric(day)<10,paste0('0',day),day))|>
    dplyr::mutate(month=monthlist[month])|>
    tidyr::unite(year,month,day,col='date',sep='-')|>
    dplyr::mutate(min=gsub('PM',':00',min))|>
    dplyr::mutate(hour=ifelse(hour=='12','00',hour))|>
    dplyr::mutate(hour=as.numeric(hour)+12)|>
    tidyr::unite(hour,min,sep=':',col='time')|>
    tidyr::unite(date,time,col=newdate,sep=' ')
  
    newtab<-tab[c(mornings,evenings,other),]
    newtab$date_sampled<-c(morntab$newdate,evetab$newdate,newtab$date_sampled[other])
    
}