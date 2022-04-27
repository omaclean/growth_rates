
library(nlstools);library(pheatmap)
new_tab_in=read.table('C:/Users/Oscar/Documents/S_for_safehaven_post20200701.tsv',header=T)
antigenic_scores=read.csv('C:/Users/Oscar/Downloads/2021-09-29_antigenic_support.csv',header=T)
antigenic_scores$mutation=gsub('X','-',antigenic_scores$mutation)
days_to_analyse=50
seqs_to_check=1000
new_tab_in$sample_date=as.Date(new_tab_in$sample_date,origin='1970-01-01')
range(new_tab$sample_date)

new_tab_in=new_tab_in[order(new_tab_in$sample_date),]

weeks_to_check=40

dat_print=matrix(ncol=6,nrow=0)
colnames(dat_print)=c('week','mutation','growth','CI1','CI2','parent_mut')
for(week_try in 0:(weeks_to_check-1)){
  date_try=max(as.numeric(new_tab_in$sample_date))-week_try*7
  new_tab=new_tab_in[as.numeric(new_tab_in$sample_date)%in%(date_try-days_to_analyse+1):
                       (date_try),]
  
  #filter out recent data if less than 5 sequences are found in those days samples (will never leave a gap of days)
  rem_dat=c()
  for(i in 0:5){
    if(length(which(new_tab$sample_date==max(new_tab$sample_date)-i))<20){
      rem_dat=c(rem_dat,max(new_tab$sample_date-i))
    }else{break}
  }
  print(as.Date(rem_dat,origin='1970-01-01'))
  
  new_tab=new_tab[!new_tab$sample_date%in%rem_dat,]
  print(c('weektry',week_try,'dim:',dim(new_tab)))
  ################################################################################################
  ################################################################################################
  ################################################################################################
  ################################################################################################
  ################################################################
  ################################
  get_muts=tail(new_tab$s,500)
  muts=as.character((unlist(sapply(get_muts,function(x) strsplit(x,';')[[1]]))))
  muts=muts[sapply(muts,function(x) substr(x,1,1)!=substr(x,nchar(x),nchar(x)))]
  muts=muts[!grepl('X',muts)]
  
  muts_tab=table(muts)
  #>1 %
  muts=names(muts_tab)[muts_tab>seqs_to_check*0.01]
  #add in those of interest if not already present
  #muts=unique(c('R346K','V213G',"A701V",muts))
  
  muts_dropout=as.character(sapply(muts,function(x) paste(substr(x,1,(nchar(x)-1)),'X',sep='')))
  names(muts_dropout)=muts
  #need to keep samples with failures to have any power due to omicron primer failure
  #new_tab=new_tab[!grepl(paste(muts_dropout,collapse='|'),new_tab$s),]
  
  #order by location in spike
  muts=muts[order(as.numeric(sapply(muts,function(x) substr(x,2,nchar(x)-1))))]
  
  muts_to_test=muts
  
  ## calculate how many of each mutations there were each day, and how many times that site failed to be called
  new_tab2=new_tab_fails=new_tab[,1:3]
  i=1
  for(i in 1:length(muts_to_test)){
    new_tab2[,muts_to_test[i]]=grepl(muts_to_test[i],new_tab$s)
    new_tab_fails[,muts_to_test[i]]=grepl(muts_dropout[i],new_tab$s)
  }
  #colnames(new_tab2)=c(colnames(new_tab2)[1:3],muts_to_test)
  
  out_tab=table(new_tab2$sample_date)
  out_tab=out_tab_fails=data.frame(date=as.Date(names(out_tab),origin='1970-01-01'),
                                   sequence_count=as.numeric(out_tab))
  for(i in 1:length(muts_to_test)){
    out_tab[,muts_to_test[i]]=table(new_tab$sample_date,
                                    grepl(muts_to_test[i],
                                          new_tab$s))[,2]
    out_tab_fails[,muts_to_test[i]]=table(new_tab$sample_date,
                                          grepl(paste(substr(muts_to_test[i],1,nchar(muts_to_test[i])-1),'X',collapse='',sep=''),
                                                new_tab$s))[,2]
  }
  
  #################################### model growth
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  for(i in 1:length(muts_to_test)){
    dat=tail(data.frame(date=out_tab$date,sequence_count=out_tab$sequence_count-out_tab_fails[,muts_to_test[i]],
                        mut=as.numeric(out_tab[,muts_to_test[i]])),days_to_analyse)
    dat$prop=dat$mut/dat$sequence_count
    
    growth_init=summary(lm(log(mut/sequence_count+1)~log(32+as.numeric(date)-max(as.numeric(dat$date))),data=tail(dat,20)))$coef[2,1]
    
    # plot(dat$date,dat$mut/dat$sequence_count)
    mod=tryCatch(nls(mut~sequence_count/(1+exp(-growth*(as.numeric(date)-constant))),data=dat,
                     start=list(growth=abs(growth_init),
                                constant=tail(as.numeric(dat$date),1)+1/tail(dat$prop,1)),
                     algorithm='port'),
                 error=function(e){return('model.fail')})
    if(mod[1]=='model.fail'){
      print(c(sum(dat$mut),muts_to_test[i]))
      next
    }
    dat_print=rbind(dat_print,c(date_try,muts_to_test[i],round(summary(mod)$coef[1,1],3),round(confint2(mod)[1,1],3),
                                round(confint2(mod)[1,2],3),c('sequence_count')))
    

  }
  
  
}

print(dat_print)





colnames(dat_print)=c('week','mutation','growth','CI1','CI2','parent_mut')
dat_print2=dat_print
library(plotrix)
dat_print=as.data.frame(dat_print)
dat_print[,1]=as.Date(as.numeric(dat_print[,1]),origin='1970-01-01')
dat_print[,3]=as.numeric(dat_print[,3])
dat_print[,4]=as.numeric(dat_print[,4])
dat_print[,5]=as.numeric(dat_print[,5])


xlim_plot=as.numeric(c(max(as.numeric(new_tab_in$sample_date))-weeks_to_check*7,
                       max(as.numeric(new_tab_in$sample_date))+2))

rows=8
par(mfrow=c(ceiling(length(unique(dat_print[,2]))/rows),rows),
    mar=c(1,1,1,0))
# par(mfrow=c(1,1),
#     mar=c(1,1,1,0))

len=length(unique(dat_print[,2]));count=0

tab_hits=table(dat_print[,2])
len=(length(tab_hits[tab_hits>2]))
par(mfrow=c(ceiling(length(tab_hits[tab_hits>2])/rows),rows),
    mar=c(2,1,1,0),new=F)
ylim_plot=as.numeric(range(c(dat_print[,4],dat_print[,5])))
ylim_plot=c(-.32,.4)
count=0

hits=names(tab_hits[tab_hits>2])
hits=hits[order(as.numeric(gsub('[A-Z]|-','',hits)))]
for(mut_hits in hits){
  count=count+1
  if(count%%rows==1){
    par(mar=c(1,2,1,0))
    if(count%in%(len-rows+1):len){
      par(mar=c(2,2,1,0))
      plotCI(dat_print[dat_print[,2]==mut_hits,1],
             dat_print[dat_print[,2]==mut_hits,3],
             liw=abs(dat_print[dat_print[,2]==mut_hits,4]-dat_print[dat_print[,2]==mut_hits,3]),
             uiw=dat_print[dat_print[,2]==mut_hits,5]-dat_print[dat_print[,2]==mut_hits,3],
             xlim=xlim_plot,ylim=ylim_plot,
             main=paste(mut_hits,'of seqs'),xlab=''
      )
      abline(h=0,lty=2)
      next
      }
    plotCI(dat_print[dat_print[,2]==mut_hits,1],
           dat_print[dat_print[,2]==mut_hits,3],
           liw=abs(dat_print[dat_print[,2]==mut_hits,4]-dat_print[dat_print[,2]==mut_hits,3]),
           uiw=dat_print[dat_print[,2]==mut_hits,5]-dat_print[dat_print[,2]==mut_hits,3],
           xlim=xlim_plot,ylim=ylim_plot,
           main=paste(mut_hits,'of seqs'),xlab='',xaxt='n'
    )
    abline(h=0,lty=2)
    next
  }else{par(mar=c(1,1,1,0))}
  
  if(count%in%(len-rows+1):len){
    par(mar=c(2,1,1,0))
    plotCI(dat_print[dat_print[,2]==mut_hits,1],
           dat_print[dat_print[,2]==mut_hits,3],
           liw=abs(dat_print[dat_print[,2]==mut_hits,4]-dat_print[dat_print[,2]==mut_hits,3]),
           uiw=dat_print[dat_print[,2]==mut_hits,5]-dat_print[dat_print[,2]==mut_hits,3],
           xlim=xlim_plot,ylim=ylim_plot,
           main=paste(mut_hits,'of seqs'),xlab='',yaxt='n'
           )
  }else{
    
    plotCI(dat_print[dat_print[,2]==mut_hits,1],
           dat_print[dat_print[,2]==mut_hits,3],
           liw=abs(dat_print[dat_print[,2]==mut_hits,4]-dat_print[dat_print[,2]==mut_hits,3]),
           uiw=dat_print[dat_print[,2]==mut_hits,5]-dat_print[dat_print[,2]==mut_hits,3],
           xlim=xlim_plot,ylim=ylim_plot,
           main=paste(mut_hits,'of seqs'),xlab='',yaxt='n',xaxt='n')
    
  }
  abline(h=0,lty=2)
  
}

#write.csv(dat_print,file='C:/Users/Oscar/Documents/mutation_growth_rate_hisotric_may2021.to.jan.2022.csv')

par(mfrow=c(1,1),
    mar=c(5,5,3,3))
mut_hits='A222V'

plotCI(dat_print[dat_print[,2]==mut_hits,1],
       dat_print[dat_print[,2]==mut_hits,3],
       liw=abs(dat_print[dat_print[,2]==mut_hits,4]-dat_print[dat_print[,2]==mut_hits,3]),
       uiw=dat_print[dat_print[,2]==mut_hits,5]-dat_print[dat_print[,2]==mut_hits,3],
       xlim=xlim_plot,ylim=ylim_plot,
       main=paste(mut_hits,'of seqs',paste(as.character(as.Date(xlim_plot,
                                                                      origin='1970-01-01')),collapse=' to ')),xlab='date',ylab='growth rate')
abline(h=0,lty=2)


library(ggplot2)
ggplot(new_tab_in,aes(x=new_tab_in$sample_date,fill=lineage))+geom_histogram()

range=as.Date(c('2021-09-01','2022-01-29'),origin='1970-01-01')
range=range[1]:range[2]
lins=c('^B.1.617.2|^AY','^AY.4.2','^BA.1','^BA.1.1','^BA.2')
names=c('basal Delta','AY.4.2','BA.1','BA.1.1','BA.2')
data=data.frame(matrix(ncol=4,nrow=0))
for(i in 1:length(range)){
  for(j in 1:length(lins)){
    data=rbind(data,c(range[i],lins[j],length(grep(lins[j],new_tab_in$lineage[
      new_tab_in$sample_date==range[i]])),names[j]))
  }
}
data=as.data.frame(data)
data_save=data
data[,3]=as.numeric(data[,3])
data[,1]=as.Date(data[,1],origin='1970-01-01')
rem=c()
for(i in 1:length(range)){
  data[data[,1]==range[i]&data[,2]==lins[1],3]=data[data[,1]==range[i]&data[,2]==lins[1],3]-
    data[data[,1]==range[i]&data[,2]==lins[2],3]
  data[data[,1]==range[i]&data[,2]==lins[3],3]=data[data[,1]==range[i]&data[,2]==lins[3],3]-
    data[data[,1]==range[i]&data[,2]==lins[4],3]
  if(sum(data[data[,1]==range[i],3])<10){print(range[i]);rem=c(rem,range[i])}
  data[data[,1]==range[i],3]=data[data[,1]==range[i],3]/sum(data[data[,1]==range[i],3])
}
head(data)
colnames(data)=c('date','lineage_regex','proportion','lineage_type')
print(as.Date(rem,origin='1970-01-01'))
data[,1]=as.Date(as.numeric(data[,1]),origin='1970-01-01')
data$lineage_type=factor(data$lineage_type,levels=names)
ggplot(data,aes(y=proportion,x=date,fill=lineage_type,color=lineage_type))+geom_bar(stat='identity',
                                                                  position='stack')+
  scale_x_date(date_labels = "%m-%Y")+theme_minimal()#+scale_fill_manual(values=(RColorBrewer::brewer.pal(length(lins),'Dark2')))+
  #scale_color_manual(values=(RColorBrewer::brewer.pal(length(lins),'Dark2')))


unique(grep('^BA.1',value=T,new_tab_in$lineage))
