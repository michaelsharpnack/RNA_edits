#this function calculates the correlation between target RNA, edits, and ADAR 1/2 RNA
#input = RNA abundance & editing frequency matrices
#output = ADAR 1/2 <-> target RNA cor matrix & ADAR1/2 <-> edits / edits <-> target RNA cor matrix
#all correlations are in spearman and then BH q-value

rna.edit.adar.cor <- function(RNA,EDIT){
  #editing frequency to ADAR1/2 correlation
  adar.edit.cor <- matrix(0,nrow=dim(EDIT)[1],ncol=4)
  rownames(adar.edit.cor) <- rownames(EDIT)
  colnames(adar.edit.cor) <- c('ADAR-EDIT Rho','ADAR-EDIT q-value','ADARB1-EDIT Rho','ADARB1-EDIT q-value')
  for(j in 1:dim(EDIT)[1]){
    if(sum(is.na(EDIT[j,]) == FALSE) > 10){
      cor.temp.1 <- cor.test(EDIT[j,][is.na(EDIT[j,]) == FALSE],RNA['ADAR',][is.na(EDIT[j,]) == FALSE],method='spearman')
      cor.temp.2 <- cor.test(EDIT[j,][is.na(EDIT[j,]) == FALSE],RNA['ADARB1',][is.na(EDIT[j,]) == FALSE],method='spearman')
      adar.edit.cor[j,1] <- cor.temp.1$estimate
      adar.edit.cor[j,2] <- cor.temp.1$p.value
      adar.edit.cor[j,3] <- cor.temp.2$estimate
      adar.edit.cor[j,4] <- cor.temp.2$p.value
    } else {
      adar.edit.cor[j,1] <- 0
      adar.edit.cor[j,2] <- 1
      adar.edit.cor[j,3] <- 0
      adar.edit.cor[j,4] <- 1
    }
  }
  adar.edit.cor[,2] <- p.adjust(adar.edit.cor[,2],method='BH')
  adar.edit.cor[,4] <- p.adjust(adar.edit.cor[,4],method='BH')
  
  #ADAR to target gene correlation
  adar.gene.cor <- matrix(0,nrow=dim(RNA)[1],ncol=4)
  rownames(adar.gene.cor) <- rownames(RNA)
  colnames(adar.gene.cor) <- c('ADAR-GENE Rho','ADAR-GENE q-value','ADARB1-GENE Rho','ADARB1-GENE q-value')
  for(j in 1:dim(RNA)[1]){
    cor.temp.1 <- cor.test(RNA[j,],RNA['ADAR',],method='spearman')
    cor.temp.2 <- cor.test(RNA[j,],RNA['ADARB1',],method='spearman')
    adar.gene.cor[j,1] <- cor.temp.1$estimate
    adar.gene.cor[j,2] <- cor.temp.1$p.value
    adar.gene.cor[j,3] <- cor.temp.2$estimate
    adar.gene.cor[j,4] <- cor.temp.2$p.value
  }
  adar.gene.cor[,2] <- p.adjust(adar.gene.cor[,2],method='BH')
  adar.gene.cor[,4] <- p.adjust(adar.gene.cor[,4],method='BH')
  
  #editing frequency to target gene correlation
  edit.rna.cor <- vector(mode='numeric',100000)
  edit.rna.cor.labels <- matrix(NA,nrow=100000,ncol=2)
  edit.rna.cor.nsamples <- vector(mode='numeric',100000)
  edit.rna.cor.mean.freq <- vector(mode='numeric',100000)
  edit.rna.cor.pvalues <- vector(mode='numeric',100000)
  k = 0
  for(i in 1:dim(EDIT)[1]){
    temp = strsplit(rownames(EDIT)[i],"\\|")[[1]]
    if(temp[3] == 'intergenic'){
      temp = strsplit(rownames(EDIT)[i],"\\(")[[1]]
      if(strsplit(temp[2],"\\,")[[1]][2] == "NONE"){
        if(length(which(rownames(RNA) == strsplit(temp[1],"\\|")[[1]][4])) == 1){
          k = k+1 
          edit.rna.cor.labels[k,1] = rownames(EDIT)[i]
          edit.rna.cor.labels[k,2] = strsplit(temp[1],"\\|")[[1]][4]
          temp10 <- try(cor.test(RNA[strsplit(temp[1],"\\|")[[1]][4],is.na(EDIT[i,]) == FALSE],
                                 EDIT[i,is.na(EDIT[i,]) == FALSE],method='spearman'))
          if(length(temp10) > 1){
            edit.rna.cor[k] = temp10$estimate
            edit.rna.cor.pvalues[k] = temp10$p.value
          }
          edit.rna.cor.nsamples[k] <- sum(is.na(EDIT[i,]) == FALSE)
          edit.rna.cor.mean.freq[k] <- mean(EDIT[i,][is.na(EDIT[i,]) == FALSE])
        }
      } else{
        if(length(which(rownames(RNA) == strsplit(temp[1],"\\|")[[1]][4])) == 1){
          k = k+1
          edit.rna.cor.labels[k,1] = rownames(EDIT)[i]
          edit.rna.cor.labels[k,2] = strsplit(temp[1],"\\|")[[1]][4]
          temp10 <- try(cor.test(RNA[strsplit(temp[1],"\\|")[[1]][4],is.na(EDIT[i,]) == FALSE],
                                 EDIT[i,is.na(EDIT[i,]) == FALSE],method='spearman'))
          if(length(temp10) > 1){
            edit.rna.cor[k] = temp10$estimate
            edit.rna.cor.pvalues[k] = temp10$p.value
          }
          edit.rna.cor.nsamples[k] <- sum(is.na(EDIT[i,]) == FALSE)
          edit.rna.cor.mean.freq[k] <- mean(EDIT[i,][is.na(EDIT[i,]) == FALSE])
          
        }
        if(length(which(rownames(RNA) == strsplit(temp[2],"\\,")[[1]][2])) == 1){
          k = k+1
          edit.rna.cor.labels[k,2] = strsplit(temp[2],"\\,")[[1]][2]
          edit.rna.cor.labels[k,1] = rownames(EDIT)[i]
          temp10 <- try(cor.test(RNA[strsplit(temp[2],"\\,")[[1]][2],is.na(EDIT[i,]) == FALSE],
                                 EDIT[i,is.na(EDIT[i,]) == FALSE],method='spearman'))
          if(length(temp10) > 1){
            edit.rna.cor[k] = temp10$estimate
            edit.rna.cor.pvalues[k] = temp10$p.value
          }
          edit.rna.cor.nsamples[k] <- sum(is.na(EDIT[i,]) == FALSE)
          edit.rna.cor.mean.freq[k] <- mean(EDIT[i,][is.na(EDIT[i,]) == FALSE])
        }
      }
    } else {
      if(length(which(rownames(RNA) == temp[4])) == 1){
        k = k+1
        edit.rna.cor.labels[k,1] = rownames(EDIT)[i]
        edit.rna.cor.labels[k,2] = temp[4]
        temp10 <- try(cor.test(RNA[temp[4],is.na(EDIT[i,]) == FALSE],
                               EDIT[i,is.na(EDIT[i,]) == FALSE],method='spearman'))
        if(length(temp10) > 1){
          edit.rna.cor[k] = temp10$estimate
          edit.rna.cor.pvalues[k] = temp10$p.value
        }
        edit.rna.cor.nsamples[k] <- sum(is.na(EDIT[i,]) == FALSE)
        edit.rna.cor.mean.freq[k] <- mean(EDIT[i,][is.na(EDIT[i,]) == FALSE])
      }
    }
    if(i %% 10000 == 0){
      print(i)
    }
  }
  edit.rna.cor <- edit.rna.cor[1:k]
  edit.rna.cor.labels <- edit.rna.cor.labels[1:k,]
  edit.rna.cor.nsamples <- edit.rna.cor.nsamples[1:k]
  edit.rna.cor.mean.freq <- edit.rna.cor.mean.freq[1:k]
  edit.rna.cor.pvalues <- edit.rna.cor.pvalues[1:k]
  edit.rna.cor.pvalues <- p.adjust(edit.rna.cor.pvalues,method='BH')
  edit.rna.cor <- cbind(edit.rna.cor.labels,edit.rna.cor,edit.rna.cor.pvalues,edit.rna.cor.mean.freq,edit.rna.cor.nsamples)
  rm(edit.rna.cor.labels,edit.rna.cor.nsamples,edit.rna.cor.mean.freq,edit.rna.cor.pvalues)
  colnames(edit.rna.cor) <- c('EDIT Name','Target Gene Symbol','EDIT-Target GENE Spearman RHO','EDIT-Target GENE q-value','EDIT mean freq','Samples compared')
  edit.rna.cor <- cbind(edit.rna.cor,adar.edit.cor[edit.rna.cor[,1],],adar.gene.cor[edit.rna.cor[,2],])
  edit.rna.cor.save <- list('edit.rna.cor' = edit.rna.cor, 'adar.edit.cor' = adar.edit.cor, 'adar.gene.cor' = adar.gene.cor)
  return(edit.rna.cor.save)
}

#edit.cor.all <- list()
#for(i in 1:4){
#  edit.cor.all[[i]] <- rna.edit.adar.cor(rna[[i]],edits[[i]])
#  print(i)
#}
