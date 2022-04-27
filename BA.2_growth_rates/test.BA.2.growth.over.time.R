
library(nlstools);library(pheatmap)
antigenic_scores=read.csv('C:/Users/Oscar/Downloads/2021-09-29_antigenic_support.csv',header=T)
antigenic_scores$mutation=gsub('X','-',antigenic_scores$mutation)

seqs_to_check=1000

dir='C:/Users/Oscar/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/oscar/genom.dat.25.3.22'
setwd(dir)

#overall window= 
weeks_to_check=4
week_start=as.Date('2021-12-01',origin='1970-01-01')
week_end=as.Date('2022-03-14',origin='1970-01-01')

time_period=list()
windows=1+ceiling(round(as.numeric(week_end-week_start))/7)-weeks_to_check
print(week_start+windows*7+weeks_to_check*7)
for(i in 1:windows){
  time_period[[i]]=c(week_start+i*7-7,week_start+(i+weeks_to_check-1)*7)
}

mut_pairs=list()
for(i in 1:length(time_period)){
  mut_pairs[[i]]=c('V213G','R346R')
}

#muts_to_plot=unique(c('R346K','V213G',"A701V","N764K",head(test_muts$muts[order(as.numeric(test_muts$grow),decreasing=T)],3)))


countries= list.files(dir,pattern='.tsv')
all_dat=matrix(ncol=6,nrow=0)
for(country in countries){
  
  upper=lower=vals=rep(0,length(time_period))
  par(mfrow=c(ceiling(sqrt(length(time_period))),ceiling(sqrt(length(time_period)))))
  
  new_tab_in=read.table(country,header=T)
  new_tab_in$sample_date=as.Date(new_tab_in$sample_date,origin='1970-01-01')
  new_tab_in=new_tab_in[order(new_tab_in$sample_date),]
  
  
  for(data_set in 1:length(mut_pairs)){
    time_window=as.Date(time_period[[data_set]],origin='1970-01-01')
    mut_pair=mut_pairs[[data_set]]
    muts_to_plot=mut_pair[1]
    muts_to_plot_rel=mut_pair[2]
    new_tab=new_tab_in[as.numeric(new_tab_in$sample_date)%in%as.numeric(time_window[1]):
                         as.numeric(time_window[2]),]
    
    if(!any(grepl(muts_to_test[i],new_tab$s))){next}
    new_tab2=new_tab[,1:3]
    i=1
    muts_to_test=mut_pair[mut_pair!='sequence_count']
    muts_dropout=c()
    for(i in 1:length(muts_to_test)){
      new_tab2[,muts_to_test[i]]=grepl(muts_to_test[i],new_tab$s)
      muts_dropout=c(muts_dropout,
                     paste(substr(muts_to_test[i],1,nchar(muts_to_test[i])-1),'X',
                           collapse='',sep=''))
    }
    names(muts_dropout)=muts_to_test
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
    
    dat=#tail(
      data.frame(date=out_tab$date,sequence_count=out_tab[,muts_to_plot_rel]-
                   out_tab_fails[,muts_to_plot[1]],
                 mut=as.numeric(out_tab[,muts_to_plot[1]]))#,days_to_analyse)
    
    
    
    #will miss sequencing errors in the containing sample- if they're numerous will be an issue, 
    #should ideally go back to the original data here as coded underneath, but it's slower to do it properly
    dat$sequence_count=sapply(1:nrow(dat),function(x) max(c(dat$mut[x],dat$sequence_count[x])))
    
    if(muts_to_plot_rel!='sequence_count'){
      dat$mut=table(new_tab$sample_date,grepl(muts_to_plot_rel,new_tab$s)&grepl(muts_to_plot[1],new_tab$s))[,2]
      
      dat$sequence_count=table(new_tab$sample_date,grepl(muts_to_plot_rel,new_tab$s)&(!grepl(muts_dropout[muts_to_plot_rel],new_tab$s)&
                                                                                        !grepl(muts_dropout[muts_to_plot[1]],new_tab$s)))[,2]
    }
    dat$prop=dat$mut/dat$sequence_count
    
    growth_init=summary(lm(log(mut/sequence_count+1)~log(32+as.numeric(date)-max(as.numeric(dat$date))),data=tail(dat,20)))$coef[2,1]
    mod='model.fail'
    trials=0
    while(mod[1]=='model.fail'|trials<10){
      mod=tryCatch(nls(mut~sequence_count/(1+exp(-growth*(as.numeric(date)-constant))),data=dat,
                       start=list(growth=abs(growth_init),
                                  constant=tail(as.numeric(dat$date),1)-7+trials+(2/(.1+tail(dat$prop,1)))),algorithm='port'),
                   error=function(e){return('model.fail')})
      
      trials=trials+1
    }
    if(mod[1]=='model.fail'){
      plot(dat$date,dat$mut/dat$sequence_count,main=paste('could not fit',hits_plot_name),
           ylab=paste('prop of ',muts_to_plot_rel,' with',muts_to_plot))
      next
    }
    print(c(muts_to_plot,round(summary(mod)$coef[1,1],3),'confint',round(confint2(mod)[1,],3)))
    antigenic_score=antigenic_scores$score[match(muts_to_plot,antigenic_scores$mutation)]
    max=dat$mut/dat$sequence_count
    max=max(max[!is.infinite(max)&!is.na(max)])
    plot(dat$date,dat$mut/dat$sequence_count,cex=1.2,ylab='',xlab='',ylim=c(0,max*1.05))
    par(new=T)
    plot(1/(1+exp(-summary(mod)$coef[1,1]*(as.numeric(dat$date)-summary(mod)$coef[2,1]))),
         pch=17,col=2,xlab='week',
         ylab=paste('prop of',muts_to_plot_rel,' with',muts_to_plot),
         ylim=c(0,max*1.05),type='l',lwd=2,
         
         main=paste(c(muts_to_plot,as.character(time_window[1]), 'to',
                      as.character(time_window[2]) ,'\n ',country,
                      'growth rate =' ,round(summary(mod)$coef[1,1],3),'CIs',
                      round(confint2(mod)[1,1],3),'to',
                      round(confint2(mod)[1,2],3)),
                    collapse=' '),xaxt='n')
    upper[data_set]=confint2(mod)[1,2]
    lower[data_set]=confint2(mod)[1,2]
    vals[data_set]=summary(mod)$coef[1,1]
    
    all_dat=rbind(all_dat,c(country,time_window,summary(mod)$coef[1,1],confint2(mod)[1,1],confint2(mod)[1,2]))
    
    legend(x='topleft',col=1:2,pch=c(1,15),legend=c('real frequency','modelled'),bty='n')
  }
  
}
write.csv(all_dat,paste('growth.rates.set.',weeks_to_check,'week_windows.csv'))

all_dat2=as.data.frame(all_dat)
colnames(all_dat2)=c('country','time1','time2','est','upper','lower')
all_dat2[,1]=sapply(all_dat2[,1],function(x) strsplit(x,'.tsv')[[1]][1])
all_dat2[,2]=as.Date(as.numeric(all_dat2[,2]),origin='1970-01-01')
all_dat2[,4]=round(as.numeric(all_dat2[,4]),3)
all_dat2[,5]=round(as.numeric(all_dat2[,5]),3)
all_dat2[,6]=round(as.numeric(all_dat2[,6]),3)

library(ggplot2)
ggplot(as.data.frame(all_dat2),aes(y=est,x=time1,fill=country,col=country))+
  geom_point(cex=3,position=position_dodge(width = 3))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=2,
                position=position_dodge(width = 3))+theme_bw()+
  xlab(paste('start date of',weeks_to_check,'week window'))+ylab('estimated daily growth rate')+
  ggtitle('BA.2 growth rate advantage vs BA.1.* (excluding BA.1.1) in 4 UK countries')+coord_cartesian(ylim=c(0.05,0.3))




##################################################
weeks_to_check=5
all_dat=read.csv(paste('growth.rates.set.',weeks_to_check,'week_windows.csv'),row.names = 1)

all_dat2=as.data.frame(all_dat)
colnames(all_dat2)=c('country','time1','time2','est','upper','lower')
all_dat2[,1]=sapply(all_dat2[,1],function(x) strsplit(x,'.tsv')[[1]][1])
all_dat2[,2]=as.Date(as.numeric(all_dat2[,2]),origin='1970-01-01')
all_dat2[,4]=round(as.numeric(all_dat2[,4]),3)
all_dat2[,5]=round(as.numeric(all_dat2[,5]),3)
all_dat2[,6]=round(as.numeric(all_dat2[,6]),3)

library(ggplot2)
ggplot(as.data.frame(all_dat2),aes(y=est,x=time1,fill=country,col=country))+
  geom_point(cex=3,position=position_dodge(width = 3))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=2,
                position=position_dodge(width = 3))+theme_bw()+
  xlab(paste('start date of',weeks_to_check,'week window'))+ylab('estimated daily growth rate')+
  ggtitle('BA.2 growth rate advantage vs BA.1.* (excluding BA.1.1) in 4 UK countries')+coord_cartesian(ylim=c(0.05,0.3))


