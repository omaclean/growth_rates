
library(nlstools);library(pheatmap)
new_tab_in=read.table('C:/Users/Oscar/Documents/S_for_safehaven_post20200701.tsv',header=T)
antigenic_scores=read.csv('C:/Users/Oscar/Downloads/2021-09-29_antigenic_support.csv',header=T)
antigenic_scores$mutation=gsub('X','-',antigenic_scores$mutation)
days_to_analyse=50
seqs_to_check=1000
new_tab_in$sample_date=as.Date(new_tab_in$sample_date,origin='1970-01-01')
range(new_tab_in$sample_date)

new_tab_in=new_tab_in[order(new_tab_in$sample_date),]

weeks_to_check=20
variants=list('alpha','delta','AY.4.2','BA.1')
time_period=list(c('2020-11-01','2021-01-30'),c('2021-04-01','2021-06-30'),
                 c('2021-07-01','2021-10-31'),c('2021-12-01','2021-12-31'))
#'A222V', 'N501Y'
mut_pairs=list(c('N501Y','sequence_count'),c('L452R','sequence_count'),c('A222V','L452R'),
               c('N501Y','sequence_count'))




#muts_to_plot=unique(c('R346K','V213G',"A701V","N764K",head(test_muts$muts[order(as.numeric(test_muts$grow),decreasing=T)],3)))

par(mfrow=c(2,2))
for(data_set in 1:length(mut_pairs)){
  
    time_window=as.Date(time_period[[data_set]],origin='1970-01-01')
    mut_pair=mut_pairs[[data_set]]
    muts_to_plot=mut_pair[1]
    muts_to_plot_rel=mut_pair[2]
    new_tab=new_tab_in[as.numeric(new_tab_in$sample_date)%in%as.numeric(time_window[1]):
                         as.numeric(time_window[2]),]
    
    new_tab2=new_tab[,1:3]
    i=1
    muts_to_test=mut_pair[mut_pair!='sequence_count']
    muts_dropout=c()
    for(i in 1:length(muts_to_test)){
      new_tab2[,muts_to_test[i]]=grepl(muts_to_test[i],new_tab$S)
      muts_dropout=c(muts_dropout,
                     paste(substr(muts_to_test[i],1,nchar(muts_to_test[i])-1),'X',collapse='',sep=''))
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
                      as.character(time_window[2]) ,'\nAntigenic_score=',
                      antigenic_score,
                      '| growth rate =' ,round(summary(mod)$coef[1,1],3),'CIs',
                      round(confint2(mod)[1,1],3),'to',
                      round(confint2(mod)[1,2],3)),
                    collapse=' '),xaxt='n')
    
    legend(x='topleft',col=1:2,pch=c(1,15),legend=c('real frequency','modelled'),bty='n')
}
