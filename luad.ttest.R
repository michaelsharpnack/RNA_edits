

ttester <- function(data,group,name){
  ttest.mat <- matrix(nrow=dim(data)[1],ncol=5)
  rownames(ttest.mat) <- rownames(data)
  for(i in 1:dim(data)[1]){
    if(length(data[i,group==1][is.na(data[i,group==1]) == FALSE]) > 5 & length(data[i,group==2][is.na(data[i,group==2]) == FALSE]) > 5){
      ttest.temp <- try(t.test(data[i,group==1][is.na(data[i,group==1]) == FALSE],data[i,group==2][is.na(data[i,group==2]) == FALSE]))
      if(length(ttest.temp) > 1){
        ttest.mat[i,1] <- ttest.temp$statistic
        ttest.mat[i,2] <- ttest.temp$estimate[[1]]
        ttest.mat[i,3] <- ttest.temp$estimate[[2]]
        ttest.mat[i,4] <- ttest.temp$p.value
      }
    }
  }
  ttest.mat[,5] <- p.adjust(ttest.mat[,4],method='BH')
  colnames(ttest.mat) <- c(paste(name,'t-statistic'),paste(name,'tumor mean'),paste(name,'normal mean'),paste(name,'p-value'),paste(name,'BH q-value'))
  return(ttest.mat)
}

tumornormal <- function(rna.tumor,rna.normal,edit.tumor,edit.normal,nsamples,name){
  #overall differential expression
  tt.edit  <- ttester(cbind(edit.tumor,edit.normal),c(rep(1,nsamples),rep(2,nsamples)),paste(name,'edit'))
  tt.rna  <- ttester(cbind(rna.tumor,rna.normal),c(rep(1,nsamples),rep(2,nsamples)),paste(name,'rna'))
  #editing high
  adar.high <- (colMeans(edit.tumor,na.rm=TRUE)-colMeans(edit.normal,na.rm=TRUE)) > 0.01
  tt.edit.adar.high <- ttester(cbind(edit.tumor[,adar.high],edit.normal[,adar.high]),c(rep(1,sum(adar.high)),rep(2,sum(adar.high))),paste(name,'edit'))
  tt.rna.adar.high <- ttester(cbind(rna.tumor[,adar.high],rna.normal[,adar.high]),c(rep(1,sum(adar.high)),rep(2,sum(adar.high))),paste(name,'edit'))
  #editing low
  adar.low <- (colMeans(edit.tumor,na.rm=TRUE)-colMeans(edit.normal,na.rm=TRUE)) <= 0.01
  tt.edit.adar.low <- ttester(cbind(edit.tumor[,adar.low],edit.normal[,adar.low]),c(rep(1,sum(adar.low)),rep(2,sum(adar.low))),paste(name,'rna'))
  tt.rna.adar.low <- ttester(cbind(rna.tumor[,adar.low],rna.normal[,adar.low]),c(rep(1,sum(adar.low)),rep(2,sum(adar.low))),paste(name,'rna'))
  
  return(list('tt.edit' = tt.edit, 'tt.rna' = tt.rna, 
              'tt.edit.adar.high' = tt.edit.adar.high, 'tt.rna.adar.high' = tt.rna.adar.high,
              'tt.edit.adar.low' = tt.edit.adar.low, 'tt.rna.adar.low' = tt.rna.adar.low
              ))
}

#<- function(ttest.save,edit.cor){
#  edit.cor.rna <- cbind(edit.cor,ttest.save$tt.edit[rownames(edit.cor),],ttest.save$tt.rna[edit.cor[,2],])
#  edit.cor.rna.filtered <- edit.cor.rna[edit.cor.rna[,19] < 0.1 & edit.cor.rna[,24] < 0.1,][rowSums(is.na(edit.cor.rna[edit.cor.rna[,19] < 0.1 & edit.cor.rna[,24] < 0.1,])) == 0,]
#}

####################################################################################################################################

#pg.ttest <- tumornormal(pg.tumor,pg.normal,pg.edit,pg.edit.normal,41,'pg')
#tcga.names <- intersect(intersect(intersect(colnames(tcga.rna),colnames(temp.rna.normal)),colnames(tcga.edit)),colnames(temp.edit.normal))
#tcga.ttest <- tumornormal(tcga.rna[,tcga.names],temp.rna.normal[,tcga.names],tcga.edit[,tcga.names],temp.edit.normal[,tcga.names],57,'tcga')
#edit.cor <- edits.sig.results[[9]][[3]]

####################################################################################################################################







