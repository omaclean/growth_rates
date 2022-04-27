
library(nlstools);library(pheatmap);library(RColorBrewer)

#new_tab_in=read.table('C:/Users/Oscar/Downloads/Scottish_20220123_date_lineage_S.tsv',header=T)
 new_tab_in=read.table('C:/Users/Oscar/Documents/S_for_safehaven_post20200701.tsv',header=T)
antigenic_scores=read.csv('C:/Users/Oscar/Downloads/2021-09-29_antigenic_support.csv',header=T)
antigenic_scores$mutation=gsub('X','-',antigenic_scores$mutation)
days_to_analyse=50
#check this many sequences for mutations to test
seqs_to_check=1000
#min frequency of mutations to test
min_frequency=0.01
new_tab_in$sample_date=as.Date(new_tab_in$sample_date,origin='1970-01-01')
range(new_tab$sample_date)

new_tab_in=new_tab_in[order(new_tab_in$sample_date),]

new_tab=new_tab_in[as.numeric(new_tab_in$sample_date)%in%(max(as.numeric(new_tab_in$sample_date))-days_to_analyse+1):
                     (as.numeric(max(new_tab_in$sample_date))),]
#filter out recent data if less than 5 sequences are found in those days samples (will never leave a gap of days)
rem_dat=c()
for(i in 0:4){
  if(length(which(new_tab$sample_date==max(new_tab$sample_date-i))<5)){
    rem_dat=c(rem_dat,max(new_tab$sample_date-i))
  }else{break}
}
new_tab=new_tab[!new_tab$sample_date%in%rem_dat,]

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################
################################
get_muts=tail(new_tab$s,seqs_to_check)
muts=as.character((unlist(sapply(get_muts,function(x) strsplit(x,';')[[1]]))))
muts=muts[sapply(muts,function(x) substr(x,1,1)!=substr(x,nchar(x),nchar(x)))]
muts=muts[!grepl('X',muts)]

muts_tab=table(muts)
#>1 %
muts=names(muts_tab)[muts_tab>seqs_to_check*min_frequency]
#add in those of interest if not already present
muts=unique(c('R346K','V213G',"A701V",muts))
#sequencing errors for sites ending in X
muts_dropout=as.character(sapply(muts,function(x) paste(substr(x,1,(nchar(x)-1)),'X',sep='')))
names(muts_dropout)=muts
#need to keep samples with failures to have any power due to omicron primer failure
#new_tab=new_tab[!grepl(paste(muts_dropout,collapse='|'),new_tab$S),]

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

# table(new_tab$sample_date,grepl('G339D',new_tab$s))
# table(new_tab$sample_date,grepl('G339X',new_tab$s))


test_muts=data.frame(muts=muts_to_test,CI1=rep(NA,length(muts_to_test)),
                                               CI2=rep(NA,length(muts_to_test)),grow=rep(NA,length(muts_to_test)))

par(mfrow=c(length(muts_to_test),1),mar=c(4,4,1,1))


#correlation of mutation frequencies per day
cor_mat=cor(out_tab[,3:ncol(out_tab)])
#correlation coexistence of mutations in the same samples 
coex_mat=cor(new_tab2[,4:ncol(new_tab2)])
#what proportion of row mutation observations are found on column mutation background (e.g. D614G will have 100% for others but low for others on itself)
coex_mat_assym=matrix(ncol=length(muts_to_test),nrow=length(muts_to_test),dat=NA)
colnames(coex_mat_assym)=row.names(coex_mat_assym)=muts_to_test
for(i in 1:length(muts_to_test)){
  for(j in 1:length(muts_to_test)){
    coex_mat_assym[i,j]=sum(new_tab2[,3+i]&new_tab2[,3+j])/sum(new_tab2[,3+i])
  #  if(i==j){coex_mat_assym[i,j]=NA}
  }
}

#plot out correlation matrices
pheatmap(coex_mat_assym,cluster_cols=F,cluster_rows=F,na_col = '#47AF47',col=
           colorRampPalette(c('#FFFFFF',rev(colorRampPalette(brewer.pal(5,'RdYlBu')
                                                    [c(1,3,5,5)],)(200)),'#000000'))(1000)
         ,main='mutation linkage- how often is row mutation found on col mutation background')

 pheatmap(cor_mat,cluster_cols=F,cluster_rows=F,na_col = '#000000',
          col=colorRampPalette(c('#FFFFFF',rev(colorRampPalette(brewer.pal(5,'RdYlBu')
                                                            [c(1,3,5,5)],)(200)),'#000000'))(1000),
          main='frequency correlation matrix')
 pheatmap(coex_mat,cluster_cols=F,cluster_rows=F,na_col = '#000000',
          col=colorRampPalette(c('#FFFFFF',rev(colorRampPalette(brewer.pal(5,'RdYlBu')
                                                            [c(1,3,5,5)],)(200)),'#000000'))(1000)
          ,main='binary per sample correlation matrix')
#########################################################################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################


coex_mat[,20:22]

link_set=list()
hit=0
for(i in 1:(length(muts_to_test)-1)){
  for(j in (i+1):length(muts_to_test)){
    #90% correlated 
    #completely arbitrary heuristic thresholds- change to what works if it doesn't
    if(cor_mat[i,j]>.95&coex_mat[i,j]>0.5){
      hit=hit+1
      link_set[[hit]]=c(muts_to_test[i],muts_to_test[j])
      #link_set[[paste(i,j)]]=c(muts_to_test[i],muts_to_test[j])
      # if(any(muts_to_test[i]%in%link_set|muts_to_test[j]%in%link_set)){
      #   link_set[muts_to_test[i]%in%link_set|muts_to_test[j]%in%link_set]=
      #     unique(c(link_set[muts_to_test[i]%in%link_set|muts_to_test[j]%in%link_set],
      #              muts_to_test[i],muts_to_test[j]))
      # }else{
      #   link_set[[paste(i,j)]]=c(muts_to_test[i],muts_to_test[j])
      # }
    }
  }
}
par(mfrow=c(length(link_set),1),mar=c(2,3,1,1))
#merge all categories
while(length(as.character(unlist(link_set)))!=length(unique(as.character(unlist(link_set))))){
  escape=F
  for(i in 1:(length(link_set)-1)){
    for(j in (i+1):length(link_set)){
     # if(any("R346K"%in%link_set[[i]])){print(i)}
      if(any(link_set[[i]]%in%link_set[[j]])){
        link_set[[i]]=unique(c(link_set[[i]],link_set[[j]]))
        #link_set=linkset[-j]
        link_set[[j]]<-NULL
        names(link_set)=1:length(link_set)
        escape=T
        break
      }
    }
    if(escape==T){break}
  }
}
#list of mutations which are linked together. Plot together to avoid redudancy
link_set
cor_mat[,'N501Y']
par(mfrow=c(length(link_set),1),mar=c(2,3,1,1))

##plot linked links sets together - double check it's working
#par(mfrow=c(1,1),mar=c(2,3,1,1))
for(i in 1:length(link_set)){
  link_set[[i]]=link_set[[i]][order(as.numeric(sapply(link_set[[i]],function(x) substr(x,2,nchar(x)-1))))]
  for(j in 1:length(link_set[[i]])){
    plot(out_tab$date,out_tab[,link_set[[i]][j]]/(out_tab$sequence_count-
                                                                    out_tab_fails[,link_set[[i]][j]])
                                                                    ,ylab='',xlab='',
         ylim=c(0,1),type='l',col=scales::alpha(j,.3),lwd=2, main=paste('linkage group',i),
         xlim=c(range(out_tab$date)[1]-10,range(out_tab$date)[2])) 
    par(new=T)
  }
  if(length(link_set[[i]])<=15){legend(x='topleft',col=1:length(link_set[[i]]),legend=link_set[[i]],bty='n',lwd=1)
  }else{
    legend(x=range(out_tab$date)[1]-11,y=1.05,col=1:length(link_set[[i]]),legend=link_set[[i]][1:15],bty='n',lwd=1)
    legend(x=range(out_tab$date)[1],y=1.05,col=1:length(link_set[[i]]),legend=link_set[[i]][16:length(link_set[[i]])],bty='n',lwd=1)
    }
     
  
  par(new=F)
}

#################################### model growth
par(mfrow=c(1,1),mar=c(4,4,1,1))
for(i in 1:length(muts_to_test)){
  dat=tail(data.frame(date=out_tab$date,sequence_count=out_tab$sequence_count-out_tab_fails[,muts_to_test[i]],
                      mut=as.numeric(out_tab[,muts_to_test[i]])),days_to_analyse)
  dat$prop=dat$mut/dat$sequence_count
  
  #32 abitrary
  growth_init=summary(lm(log(mut/sequence_count+1)~log(32+as.numeric(date)-max(as.numeric(dat$date))),data=tail(dat,20)))$coef[2,1]
  
 # plot(dat$date,dat$mut/dat$sequence_count)
  mod=tryCatch(nls(mut~sequence_count/(1+exp(-growth*(as.numeric(date)-constant))),data=dat,
                   start=list(growth=abs(growth_init),
                              constant=tail(as.numeric(dat$date),1)+1/tail(dat$prop,1)),
                   algorithm='port'),
               error=function(e){return('model.fail')})
  # while(mod[1]=='model.fail'&trials<20){
  #   mod=tryCatch(nls(mut~sequence_count/(1+exp(-growth*(as.numeric(date)-constant))),data=dat,
  #                    start=list(growth=abs(growth_init),
  #                               constant=tail(as.numeric(dat$date),1)-10+(2/(.1+tail(dat$prop,1)))),algorithm='port'),
  #                error=function(e){return('model.fail')})
  #   
  #   trials=trials+1
  # }
  
  if(mod[1]=='model.fail'){
    print(c(sum(dat$mut),muts_to_test[i]))
    next
  }
  test_muts$CI1[i]=round(confint2(mod)[1,1],3)
  test_muts$CI2[i]=round(confint2(mod)[1,2],3)
  test_muts$grow[i]=round(summary(mod)$coef[1,1],3)
}

hits=test_muts$muts[!is.na(test_muts$CI1)&test_muts$CI1>0]
for(link in 1:length(link_set)){
  if(length(which(hits%in%link_set[[link]]))>1){
    hits=hits[-sample(which(hits%in%link_set[[link]]),
                      length(which(hits%in%link_set[[link]]))-1)]
  }
}



hits_plot_name=muts_to_plot=hits=unique(c('R346K','V213G','N501Y',"A701V","N764K",hits))

antigenic_score=antigenic_scores$score[match(muts_to_plot,antigenic_scores$mutation)]
for(link in 1:length(link_set)){
  
  if(any(muts_to_plot%in%link_set[[link]])){
    hits_plot_name[which(hits%in%link_set[[link]])]=
      paste(as.character(c(hits_plot_name[which(hits%in%link_set[[link]])],
                                                          'with',paste(link_set[[link]],collapse='&'))),collapse='')
    
    antigenic_score[which(hits%in%link_set[[link]])]=paste(antigenic_score[which(hits%in%link_set[[link]])],
                                                           sum(as.numeric(antigenic_scores$score[match(link_set[[link]],antigenic_scores$mutation)][is.finite(
                                                             as.numeric(antigenic_scores$score[match(link_set[[link]],antigenic_scores$mutation)]))]))
                                                             ,sep='&sum=')
    
  }
}





#plot
###################



#muts_to_plot=unique(c('R346K','V213G',"A701V","N764K",head(test_muts$muts[order(as.numeric(test_muts$grow),decreasing=T)],3)))
par(mfrow=c(ceiling(length(muts_to_plot)/2),2),mar=c(2,4,2.5,1))
for(i in 1:length(muts_to_plot)){
  dat=tail(data.frame(date=out_tab$date,sequence_count=out_tab$sequence_count-out_tab_fails[,muts_to_plot[i]],
                      mut=as.numeric(out_tab[,muts_to_plot[i]])),days_to_analyse)
  dat$prop=dat$mut/dat$sequence_count
  
  growth_init=summary(lm(log(mut/sequence_count+1)~log(32+as.numeric(date)-max(as.numeric(dat$date))),data=tail(dat,20)))$coef[2,1]
  mod='model.fail'
  trials=0
  while(mod[1]=='model.fail'&trials<20){
    mod=tryCatch(nls(mut~sequence_count/(1+exp(-growth*(as.numeric(date)-constant))),data=dat,
                     start=list(growth=abs(growth_init),
                                constant=tail(as.numeric(dat$date),1)-10+(2/(.1+tail(dat$prop,1)))),algorithm='port'),
                 error=function(e){return('model.fail')})
    
    trials=trials+1
  }
 # mod=gam(mut~s(sequence_count),data=dat)
  #plot(mod)

  if(mod[1]=='model.fail'){
    plot(dat$date,dat$mut/dat$sequence_count,main=paste('could not fit',hits_plot_name[i]),
         ylab=paste('prop of non-errored sequences with',muts_to_plot[i]))
    next
  }
  print(c(muts_to_plot[i],round(summary(mod)$coef[1,1],3),'confint',round(confint2(mod)[1,],3)))
  
  max=dat$mut/dat$sequence_count
  max=max(max[!is.infinite(max)&!is.na(max)])
  plot(dat$date,dat$mut/dat$sequence_count,cex=1.2,ylab='',xlab='',ylim=c(0,max*1.05))
  par(new=T)
  plot(1/(1+exp(-summary(mod)$coef[1,1]*(as.numeric(dat$date)-summary(mod)$coef[2,1]))),
       pch=17,col=2,xlab='week',ylab=paste('prop of sequences with',muts_to_plot[i]),ylim=c(0,max*1.05),type='l',lwd=2,
       main=paste(c(hits_plot_name[i],'\nAntigenic_score=',
                    antigenic_score[i],
                    '| growth rate =' ,round(summary(mod)$coef[1,1],3),'CIs',round(confint2(mod)[1,1],3),'to',
                    round(confint2(mod)[1,2],3)),
                  collapse=' '),xaxt='n')
  
  
  par(new=T)
  power=function(x,y){return(x^y)}
  dat$date_num=as.numeric(dat$date)-as.numeric(dat$date[1])+1
  mod2=(lm(abs(mut/sequence_count)~date_num+power(date_num,2)+power(date_num,3),data=dat))
  summary(mod2)
  plot(dat$date,(summary(mod2)$coefficients[1,1]+
                   dat$date_num*summary(mod2)$coefficients[2,1]+
                   power(dat$date_num,2)*summary(mod2)$coefficients[3,1]+
                   power(dat$date_num,3)*summary(mod2)$coefficients[4,1]),ylim=c(0,max*1.05),type='l',col=3
       ,xaxt='n',yaxt='n',ylab='',xlab='')
  
  legend(x='topleft',col=1:3,pch=c(1,15,15),legend=c('real frequency','modelled','cubic spline fit'),bty='n')
  
}


################# daughter mutation


#recreates coex_mat_assym but more query-able
#  more than 97% of query mutation linked to other mutation, 
#    but query mutation present for <90% of instances of other mutation 
muts_to_plot_rel=rep('sequence_count',length(muts_to_plot))
for(i in 1:length(muts_to_plot)){
  other_muts=colnames(new_tab2)[4:ncol(new_tab2)]
  other_muts=other_muts[other_muts!=muts_to_plot[i]&other_muts!='D614G']
  coex=sum_other=rep(0,length(other_muts))
  names(coex)=names(sum_other)=other_muts
  for(j in 1:length(other_muts)){
    coex[j]=sum(new_tab2[,other_muts[j]]&new_tab2[,muts_to_plot[i]])
    sum_other[j]=sum(new_tab2[,other_muts[j]])
  }
  if(any(coex>0.97*sum(new_tab2[,muts_to_plot[i]])&
         sum(new_tab2[,muts_to_plot[i]])<0.9*sum_other)){
    muts_to_plot_rel[i]=names(coex)[which(coex>0.97*sum(new_tab2[,muts_to_plot[i]])&
                                          sum(new_tab2[,muts_to_plot[i]])<0.9*sum_other)[1]]
  }
}



#muts_to_plot=unique(c('R346K','V213G',"A701V","N764K",head(test_muts$muts[order(as.numeric(test_muts$grow),decreasing=T)],3)))
par(mfrow=c(ceiling(length(muts_to_plot)/2),2),mar=c(2,4,2.5,1))
for(i in 1:length(muts_to_plot)){
  dat=tail(data.frame(date=out_tab$date,sequence_count=out_tab[,muts_to_plot_rel[i]]-out_tab_fails[,muts_to_plot[i]],
                      mut=as.numeric(out_tab[,muts_to_plot[i]]),N501Y=out_tab$N501Y),days_to_analyse)
  
  
  
  #will miss sequencing errors in the containing sample- if they're numerous will be an issue, 
  #should ideally go back to the original data here as coded underneath, but it's slower to do it properly
  dat$sequence_count=sapply(1:nrow(dat),function(x) max(c(dat$mut[x],dat$sequence_count[x])))
  
  if(muts_to_plot_rel[i]!='sequence_count'){
    dat$mut=table(new_tab$sample_date,grepl(muts_to_plot_rel[i],new_tab$s)&grepl(muts_to_plot[i],new_tab$s))[,2]

    dat$sequence_count=table(new_tab$sample_date,grepl(muts_to_plot_rel[i],new_tab$s)&(!grepl(muts_dropout[muts_to_plot_rel[i]],new_tab$s)&
                                                     !grepl(muts_dropout[muts_to_plot[i]],new_tab$s)))[,2]
  }
  dat$prop=dat$mut/dat$sequence_count
  
  growth_init=summary(lm(log(mut/sequence_count+1)~log(32+as.numeric(date)-max(as.numeric(dat$date))),data=tail(dat,20)))$coef[2,1]
  mod='model.fail'
  trials=0
  while(mod[1]=='model.fail'&trials<20){
    mod=tryCatch(nls(mut~sequence_count/(1+exp(-growth*(as.numeric(date)-constant))),data=dat,
                     start=list(growth=abs(growth_init),
                                constant=tail(as.numeric(dat$date),1)-7+trials+(2/(.1+tail(dat$prop,1)))),algorithm='port'),
                 error=function(e){return('model.fail')})
    
    trials=trials+1
  }
  if(mod[1]=='model.fail'){
    plot(dat$date,dat$mut/dat$sequence_count,main=paste('could not fit',hits_plot_name[i]),
         ylab=paste('prop of ',muts_to_plot_rel[i],' with',muts_to_plot[i]))
    next
  }
  print(c(muts_to_plot[i],round(summary(mod)$coef[1,1],3),'confint',round(confint2(mod)[1,],3)))
  
  max=dat$mut/dat$sequence_count
  max=max(max[!is.infinite(max)&!is.na(max)])
  plot(dat$date,dat$mut/dat$sequence_count,cex=1.2,ylab='',xlab='',ylim=c(0,max*1.05))
  par(new=T)
  plot(1/(1+exp(-summary(mod)$coef[1,1]*(as.numeric(dat$date)-summary(mod)$coef[2,1]))),
       pch=17,col=2,xlab='week',
       ylab=paste('prop of',muts_to_plot_rel[i],' with',muts_to_plot[i]),
       ylim=c(0,max*1.05),type='l',lwd=2,
       
       main=paste(c(hits_plot_name[i],'\nAntigenic_score=',
                    antigenic_score[i],
                    '| growth rate =' ,round(summary(mod)$coef[1,1],3),'CIs',
                    round(confint2(mod)[1,1],3),'to',
                    round(confint2(mod)[1,2],3)),
                  collapse=' '),xaxt='n')
  
  legend(x='topleft',col=1:2,pch=c(1,15),legend=c('real frequency','modelled'),bty='n')
}

# table of mutation growth rates as a proportion of sequences per day.
kable(test_muts)
