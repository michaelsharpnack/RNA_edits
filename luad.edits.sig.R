#This file produces a list of genes / edits / fasta segments to search for RBP motifs / miRNA segments
#filename: edits.sig.R

####################################################################################################################################

#Main function
edits.sig <- function(edit.rna.cor.filtered,qvalue,edit.type = NULL){
  #filter by EDIT to Target correlation
  edit.rna.cor.filtered <- edit.rna.cor.filtered[as.numeric(edit.rna.cor.filtered[,4]) < qvalue,]
  #filter by EDIT mean frequency

  #filter by # Samples with >10 reads

  #filter by EDIT to ADAR1/2 RNA cor

  #filter by Target to ADAR1/2 RNA cor

  #filter by location in GENE of EDIT

  edits.annot <- matrix(NA,nrow=dim(edit.rna.cor.filtered)[1],ncol=2)
  names.temp <- strsplit(edit.rna.cor.filtered[,1],'\\|')
  for(i in 1:dim(edit.rna.cor.filtered)[1]){
    edits.annot[i,1] <- names.temp[[i]][4]
    edits.annot[i,2] <- names.temp[[i]][3]
  }
  if(is.null(edit.type) == FALSE){
    edit.rna.cor.filtered <- edit.rna.cor.filtered[which(edits.annot[,2] == edit.type),]
  }

  #Compile by genes
  sig.genes <- list()
  genes.info <- matrix(0,nrow=length(unique(edit.rna.cor.filtered[,2])),ncol=3)
  colnames(genes.info) <- c('Sig EDITs per Target','Sig edits positive cor','Sig edits negative cor')
  rownames(genes.info) <- unique(edit.rna.cor.filtered[,2])
  for(i in 1:length(unique(edit.rna.cor.filtered[,2]))){
    temp <- edit.rna.cor.filtered[which(edit.rna.cor.filtered[,2] == unique(edit.rna.cor.filtered[,2])[i]),]
    if(is.null(dim(temp)) == FALSE){
      sig.genes[[i]] <- temp[order(temp[,1]),]
      genes.info[i,1] <- dim(sig.genes[[i]])[1]
      genes.info[i,2] <- sum(sig.genes[[i]][,3] > 0)
      genes.info[i,3] <- sum(sig.genes[[i]][,3] < 0)
    } else {
      sig.genes[[i]] <- temp
      genes.info[i,1] <- 1
      if(sig.genes[[i]][3] > 0){
        genes.info[i,2] <- 1
      } else {
        genes.info[i,3] <- 1
      }
    }
  }
  edits.annot <- edits.annot[order(edit.rna.cor.filtered[,1]),]
  edit.rna.cor.filtered <- edit.rna.cor.filtered[order(edit.rna.cor.filtered[,1]),]
  edits.annot <- edits.annot[rowSums(is.na(edit.rna.cor.filtered)) == 0,]
  edit.rna.cor.filtered <- edit.rna.cor.filtered[rowSums(is.na(edit.rna.cor.filtered)) == 0,]
  return(list(edits.annot,genes.info,edit.rna.cor.filtered[order(edit.rna.cor.filtered[,1]),]))
}

#initialize input/output
#edits.sig.results <- list()
#edit.rna.cor.filtered <- edit.cor.all[[4]]$edit.rna.cor
#edits.annot <- matrix(NA,nrow=dim(edit.rna.cor.filtered)[1],ncol=2)
#names.temp <- strsplit(edit.rna.cor.filtered[,1],'\\|')
#for(i in 1:dim(edit.rna.cor.filtered)[1]){
#  edits.annot[i,1] <- names.temp[[i]][4]
#  edits.annot[i,2] <- names.temp[[i]][3]
#}
#edits.sig.ngenes <- matrix(0,nrow=length(unique(edits.annot[,2])),ncol=2)
#rownames(edits.sig.ngenes) <- sort(unique(edits.annot[,2]))
#colnames(edits.sig.ngenes) <- c('# genes with positive cor','# genes with negative cor')

#edit.type <- NULL

#run function
#for(j in 1:length(unique(edits.annot[,2]))){
#  edits.sig.results[[j]] <- edits.sig(edit.rna.cor.filtered,0.1,sort(unique(edits.annot[,2]))[j])
#  edits.sig.ngenes[j,1] <- sum(edits.sig.results[[j]][[2]][,2] > 0)
#  edits.sig.ngenes[j,2] <- sum(edits.sig.results[[j]][[2]][,3] > 0)
#}

##################################################################################################################################

#convert frequently edited regions into bed files
bed.creator <- function(genes,edit.rna.cor.filtered,thresh,edits.annot){
  bed.save.temp <- matrix(NA,nrow=length(genes)*15,ncol=6)
  k = 0
  for(i in 1:length(genes)){
    temp52 <- which(edit.rna.cor.filtered[,2] == genes[i])
    if(length(temp52) >= thresh){
      u = 1
      v = 1
      #bed.save.temp.3 <- vector(mode='numeric',length(temp52))
      #bed.save.temp.4 <- vector(mode='numeric',length(temp52)-1)
    
      for(j in 1:length(temp52)){
        u = u+1
        if(j < length(temp52)){
          bed.save.temp.2 <- as.numeric(strsplit(edit.rna.cor.filtered[,1][temp52],'\\|')[[j+1]][2])-as.numeric(strsplit(edit.rna.cor.filtered[,1][temp52],'\\|')[[j]][2])
        }
        #bed.save.temp.3[j+1] <- bed.save.temp.2+bed.save.temp.3[j]
        #bed.save.temp.4[j] <- bed.save.temp.2
        if(bed.save.temp.2 > 100 | j == length(temp52)){
          if(u > 5){
            k = k+1
            bed.save.temp[k,3] <- as.numeric(strsplit(edit.rna.cor.filtered[,1][temp52],'\\|')[[v]][2])-20
            bed.save.temp[k,c(2,5)] <- strsplit(edit.rna.cor.filtered[,1][temp52],'\\|')[[1]][c(1,5)]
            bed.save.temp[k,1] <- strsplit(edit.rna.cor.filtered[,1][temp52],'\\|')[[1]][4]
            bed.save.temp[k,4] <- as.numeric(strsplit(edit.rna.cor.filtered[,1][temp52],'\\|')[[j]][2])+20
            bed.save.temp[k,6] <- edits.annot[temp52[1],2]
          }
          v = j+1
          u = 1
        }
      }
    }
  }
  bed.save.temp <- bed.save.temp[1:k,]
  return(bed.save.temp)
}

#bed.save.tcga <- bed.creator(unique(edits.sig.results[[9]][[3]][,2]),edits.sig.results[[9]][[3]],5)


##################################################################################################################################

#visualize results
#View(edits.sig.results[[2]][order(edits.sig.results[[2]][,2],decreasing=TRUE),])
#View(edits.sig.results[[2]][order(edits.sig.results[[2]][,3],decreasing=TRUE),])

#barplot(edits.sig.results[[2]])





