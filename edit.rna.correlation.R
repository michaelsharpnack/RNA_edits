#this file loads in the editing and RNA expression data and calculates the association between editing and expression

cancers <- c('BLCA','BRCA','CESC','CRCC','GBM','HNSC','KICH','KIRC','KIRP','LGG','LIHC','LUAD','LUSC','PRAD','STAD','THCA','UCEC')

library(data.table)

#first load in editing files of choice
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/tcga_edits/')
files.edit <- dir()
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/tcga.rnaseq/')
files.rna <- dir()
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query')

temp.rna.geneid <- read.csv(paste('tcga.rnaseq/',files.rna[1],sep=""),stringsAsFactors=FALSE,sep='\t')[,1]
temp.names <- strsplit(temp.rna.geneid,"\\|")
temp.names <- temp.names[-c(1:30)]
temp.names.2 <- vector(mode='numeric',length(temp.names))
for(i in 1:length(temp.names)){
  temp.names.2[i] <- temp.names[[i]][1]
}
rna.geneid <- temp.names.2
rm(temp.names,temp.names.2,temp.rna.geneid)

edits.union.cor <- list()
for(j in 1:length(cancers)){
  
  print(paste("Reading in editing file for ",cancers[j]))
  temp.edit <- fread(paste('tcga_edits/',files.edit[j],sep=""))
  temp.edit.rownames <- temp.edit[[1]]
  temp.edit <- as.matrix(temp.edit[,-1])
  temp.edit <- temp.edit[,-dim(temp.edit)[2]]
  class(temp.edit) <- 'numeric'
  rownames(temp.edit) <- temp.edit.rownames
  colnames(temp.edit) <- substr(colnames(temp.edit),8+nchar(cancers[j]),23)
  
  print(paste("Reading in RNAseq file for ",cancers[j]))
  temp.rna <- fread(paste('tcga.rnaseq/',files.rna[j],sep=""))
  temp.rna <- as.matrix(temp.rna[,-1])
  temp.rna <- temp.rna[-c(1:30),]
  class(temp.rna) <- 'numeric'
  rownames(temp.rna) <- rna.geneid
  colnames(temp.rna) <- substr(colnames(temp.rna),1,12)
  
  temp.rna <- temp.rna[,intersect(colnames(temp.rna),colnames(temp.edit))]
  temp.edit <- temp.edit[,intersect(colnames(temp.rna),colnames(temp.edit))]

  edit.rna.cor <- vector(mode='numeric',100000)
  edit.rna.cor.labels <- matrix(NA,nrow=100000,ncol=2)
  edit.rna.cor.nsamples <- vector(mode='numeric',100000)
  edit.rna.cor.mean.freq <- vector(mode='numeric',100000)
  edit.rna.cor.pvalues <- vector(mode='numeric',100000)
  k = 0
  for(i in 1:dim(temp.edit)[1]){
    temp = strsplit(rownames(temp.edit)[i],"\\|")[[1]]
    if(temp[3] == 'intergenic'){
      temp = strsplit(rownames(temp.edit)[i],"\\(")[[1]]
      if(strsplit(temp[2],"\\,")[[1]][2] == "NONE"){
        if(length(which(rownames(temp.rna) == strsplit(temp[1],"\\|")[[1]][4])) == 1){
          k = k+1 
          edit.rna.cor.labels[k,1] = rownames(temp.edit)[i]
          edit.rna.cor.labels[k,2] = strsplit(temp[1],"\\|")[[1]][4]
          temp10 <- try(cor.test(temp.rna[strsplit(temp[1],"\\|")[[1]][4],is.na(temp.edit[i,]) == FALSE],
                        temp.edit[i,is.na(temp.edit[i,]) == FALSE],method='spearman'))
          if(length(temp10) > 1){
            edit.rna.cor[k] = temp10$estimate
            edit.rna.cor.pvalues[k] = temp10$p.value
          }
          edit.rna.cor.nsamples[k] <- sum(is.na(temp.edit[i,]) == FALSE)
          edit.rna.cor.mean.freq[k] <- mean(temp.edit[i,][is.na(temp.edit[i,]) == FALSE])
        }
      } else{
        if(length(which(rownames(temp.rna) == strsplit(temp[1],"\\|")[[1]][4])) == 1){
          k = k+1
          edit.rna.cor.labels[k,1] = rownames(temp.edit)[i]
          edit.rna.cor.labels[k,2] = strsplit(temp[1],"\\|")[[1]][4]
          temp10 <- try(cor.test(temp.rna[strsplit(temp[1],"\\|")[[1]][4],is.na(temp.edit[i,]) == FALSE],
                        temp.edit[i,is.na(temp.edit[i,]) == FALSE],method='spearman'))
          if(length(temp10) > 1){
            edit.rna.cor[k] = temp10$estimate
            edit.rna.cor.pvalues[k] = temp10$p.value
          }
          edit.rna.cor.nsamples[k] <- sum(is.na(temp.edit[i,]) == FALSE)
          edit.rna.cor.mean.freq[k] <- mean(temp.edit[i,][is.na(temp.edit[i,]) == FALSE])
          
        }
        if(length(which(rownames(temp.rna) == strsplit(temp[2],"\\,")[[1]][2])) == 1){
          k = k+1
          edit.rna.cor.labels[k,2] = strsplit(temp[2],"\\,")[[1]][2]
          edit.rna.cor.labels[k,1] = rownames(temp.edit)[i]
          temp10 <- try(cor.test(temp.rna[strsplit(temp[2],"\\,")[[1]][2],is.na(temp.edit[i,]) == FALSE],
                        temp.edit[i,is.na(temp.edit[i,]) == FALSE],method='spearman'))
          if(length(temp10) > 1){
            edit.rna.cor[k] = temp10$estimate
            edit.rna.cor.pvalues[k] = temp10$p.value
          }
          edit.rna.cor.nsamples[k] <- sum(is.na(temp.edit[i,]) == FALSE)
          edit.rna.cor.mean.freq[k] <- mean(temp.edit[i,][is.na(temp.edit[i,]) == FALSE])
          
        }
      }
      
    } else {
      if(length(which(rownames(temp.rna) == temp[4])) == 1){
        k = k+1
        edit.rna.cor.labels[k,1] = rownames(temp.edit)[i]
        edit.rna.cor.labels[k,2] = temp[4]
        temp10 <- try(cor.test(temp.rna[temp[4],is.na(temp.edit[i,]) == FALSE],
                    temp.edit[i,is.na(temp.edit[i,]) == FALSE],method='spearman'))
        if(length(temp10) > 1){
          edit.rna.cor[k] = temp10$estimate
          edit.rna.cor.pvalues[k] = temp10$p.value
        }
        edit.rna.cor.nsamples[k] <- sum(is.na(temp.edit[i,]) == FALSE)
        edit.rna.cor.mean.freq[k] <- mean(temp.edit[i,][is.na(temp.edit[i,]) == FALSE])
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
  
  edits.union.cor[[j]] <- cbind(edit.rna.cor.labels,edit.rna.cor,edit.rna.cor.nsamples,edit.rna.cor.mean.freq,edit.rna.cor.pvalues)
  
}

edits.union <- vector(mode='character',0)
for(i in 1:length(edits.union.cor)){
  edits.union <- union(edits.union,edits.union.cor[[i]][,1])
}
edits.union <- sort(edits.union)
edits.mat <- matrix(NA,nrow=length(edits.union),length(edits.union.cor))
rownames(edits.mat) <- edits.union
colnames(edits.mat) <- cancers
edits.nsamples <- matrix(NA,nrow=length(edits.union),length(edits.union.cor))
rownames(edits.nsamples) <- edits.union
colnames(edits.nsamples) <- cancers
edits.mean.freq <- matrix(NA,nrow=length(edits.union),length(edits.union.cor))
rownames(edits.mean.freq) <- edits.union
colnames(edits.mean.freq) <- cancers
edits.pvalues <- matrix(NA,nrow=length(edits.union),length(edits.union.cor))
rownames(edits.pvalues) <- edits.union
colnames(edits.pvalues) <- cancers
edits.union.genes <- edits.union
for(i in 1:length(edits.union)){
  edits.union.genes[i] <- strsplit(edits.union[i],'\\|')[[1]][4]
  if(i %% 10000 == 0){print(i)}
}

for(i in 1:length(edits.union.cor)){
  edits.mat[edits.union.cor[[i]][,1],i] <- as.numeric(edits.union.cor[[i]][,3])
  edits.nsamples[edits.union.cor[[i]][,1],i] <- as.numeric(edits.union.cor[[i]][,4])
  edits.mean.freq[edits.union.cor[[i]][,1],i] <- as.numeric(edits.union.cor[[i]][,5])
  edits.pvalues[edits.union.cor[[i]][,1],i] <- as.numeric(edits.union.cor[[i]][,6])
  edits.pvalues[,i] <- p.adjust(edits.pvalues[,i],method='BH')
}

#assignment for Kun due on Monday. Analyse these genes.

edits.mat <- edits.mat[,-c(4,7,17)]
edits.nsamples <- edits.nsamples[,-c(4,7,17)]
edits.mean.freq <- edits.mean.freq[,-c(4,7,17)]
edits.pvalues <- edits.pvalues[,-c(4,7,17)]
View(sort(rowSums(edits.pvalues < 0.1,na.rm=TRUE),decreasing=TRUE))
View(edits.mat[order(rowSums(edits.pvalues < 0.1,na.rm=TRUE),decreasing=TRUE),])
hist(rowSums(edits.pvalues < 0.1,na.rm=TRUE),14,col='blue')



#enrichment in chromosomes/regions?
edit.temp.5 <- rowSums(edits.pvalues < 0.1,na.rm=TRUE)
edit.temp.5split <- strsplit(names(edit.temp.5),'\\|')
edit.temp.5chr <- vector(mode='character',length(edit.temp.5))
edit.temp.5bp <- edit.temp.5chr
edit.temp.5type <- edit.temp.5chr
edit.temp.5gene <- edit.temp.5chr
for(i in 1:length(edit.temp.5)){
  edit.temp.5chr[i] <- edit.temp.5split[[i]][1]
  edit.temp.5bp[i] <- as.numeric(edit.temp.5split[[i]][2])
  edit.temp.5type[i] <- edit.temp.5split[[i]][3]
  edit.temp.5gene[i] <- edit.temp.5split[[i]][4]
}
edit.temp.5chr <- gsub('chr','',edit.temp.5chr)
edit.temp.5chr[edit.temp.5chr == "X"] = 23
edit.temp.5chr[edit.temp.5chr == 'Y'] = 24
edit.temp.5chr[edit.temp.5chr == 'M'] = 25
edit.temp.5chr <- as.numeric(edit.temp.5chr)
edit.temp.5bp <- as.numeric(edit.temp.5bp)
x <- data.frame('CHR' = edit.temp.5chr, 'BP' = edit.temp.5bp, 'P' = edit.temp.5)
manhattan(x,logp=FALSE,main='Number of cancers with significant RNA-Editing correlation')

#enrichment in editing site specificity?
types <- unique(edit.temp.5type)
expect.pct <- matrix(0,nrow=length(types)+1,ncol=14)
types.n.enrich.all <- expect.pct[-1,]
for(u in 1:14){
  types.n <- vector(mode='numeric',length(types))
  names(types.n) <- types
  types.n.enrich <- types.n
  for(i in 1:length(types)){
    types.n[i] <- sum(edit.temp.5type == types[i])
    types.n.enrich[i] <- sum(edit.temp.5type == types[i] & rowSums(edits.pvalues < 0.1,na.rm=TRUE) > u-1)
  }
  types.n.expect <- types.n*sum(types.n.enrich)/sum(types.n)
  types.n.enrich.all[,u] <- types.n.enrich
  expect.pct[,u] <- c(types.n.expect[1]/types.n[1],types.n.enrich/types.n)
  #View(cbind(types.n,types.n.enrich,types.n.enrich/types.n,types.n.expect,types.n.expect/types.n)[order(names(types.n)),])
}
matplot(t(expect.pct),type='l',lwd=3,col=c(1,2,5,4,4,3,3,3,4,5,2),lty=c(1,1,1,1,2,2,1,3,3,2,2))
legend(10,0.6,c('Expected',types),cex=1,col=c(1,2,5,4,4,3,3,3,4,5,2),lty=c(1,1,1,1,2,2,1,3,3,2,2),lwd=3)

matplot(t(expect.pct),type='l',lwd=3,col=c(1,2,5,4,4,3,3,3,4,5,2),lty=c(1,1,1,1,2,2,1,3,3,2,2),ylim=c(0,0.025))
legend(12,0.025,c('Expected',types),cex=0.7,col=c(1,2,5,4,4,3,3,3,4,5,2),lty=c(1,1,1,1,2,2,1,3,3,2,2),lwd=3)

#################################################################################################################################################

#this segment produces a genewise values that you can filter on to find the most important edited genes.
#4 filtering criteria
#1. ADAR1/2 correlation with target gene
#2. ADAR1/2 correlation with edits within target gene
#3. frequency of edits within target gene
#4. edits correlation with target gene RNA


#Creates a matrix with rows edits and columns cancers. Each value is the mean editing frequency across patients within a cancer
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/tcga_edits/')
files.edit <- dir()
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/tcga.rnaseq/')
files.rna <- dir()
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query')
files.edit <- files.edit[-c(4,7,17)]
files.rna <- files.rna[-c(4,7,17)]
edits.freq <- matrix(0,nrow=dim(edits.mat)[1],ncol=14)
rownames(edits.freq) <- rownames(edits.mat)
for(j in 1:length(cancers)){
  
  print(paste("Reading in editing file for ",cancers[j]))
  temp.edit <- fread(paste('tcga_edits/',files.edit[j],sep=""))
  temp.edit.rownames <- temp.edit[[1]]
  temp.edit <- as.matrix(temp.edit[,-1])
  temp.edit <- temp.edit[,-dim(temp.edit)[2]]
  class(temp.edit) <- 'numeric'
  rownames(temp.edit) <- temp.edit.rownames
  colnames(temp.edit) <- substr(colnames(temp.edit),8+nchar(cancers[j]),23)
  
  print(paste("Reading in RNAseq file for ",cancers[j]))
  temp.rna <- fread(paste('tcga.rnaseq/',files.rna[j],sep=""))
  temp.rna <- as.matrix(temp.rna[,-1])
  temp.rna <- temp.rna[-c(1:30),]
  class(temp.rna) <- 'numeric'
  rownames(temp.rna) <- rna.geneid
  colnames(temp.rna) <- substr(colnames(temp.rna),1,12)
  
  temp.rna <- temp.rna[,intersect(colnames(temp.rna),colnames(temp.edit))]
  temp.edit <- temp.edit[,intersect(colnames(temp.rna),colnames(temp.edit))]
  
  edits.freq[intersect(rownames(temp.edit),rownames(edits.freq)),j] <- rowMeans(temp.edit[intersect(rownames(temp.edit),rownames(edits.freq)),],na.rm=TRUE)
}

#create ADAR-target gene RNA matrix, where each column is a cancer and each row is a gene
adar.rna.cor.mat <- matrix(0,nrow=length(unique(edit.temp.5gene.temp)),ncol=14)
rownames(adar.rna.cor.mat) <- unique(edit.temp.5gene.temp)
adar2.rna.cor.mat <- matrix(0,nrow=length(unique(edit.temp.5gene.temp)),ncol=14)
rownames(adar2.rna.cor.mat) <- unique(edit.temp.5gene.temp)
adar.rna.cor.mat.p <- matrix(0,nrow=length(unique(edit.temp.5gene.temp)),ncol=14)
rownames(adar.rna.cor.mat.p) <- unique(edit.temp.5gene.temp)
adar2.rna.cor.mat.p <- matrix(0,nrow=length(unique(edit.temp.5gene.temp)),ncol=14)
rownames(adar2.rna.cor.mat.p) <- unique(edit.temp.5gene.temp)
adar12.cor <- vector(mode='numeric',14)
names(adar12.cor) <- cancers
for(i in 1:14){
  adar.rna.cor.mat[,i] <- adar.rna.cor[[i]][unique(edit.temp.5gene.temp),1]
  adar2.rna.cor.mat[,i] <- adar2.rna.cor[[i]][unique(edit.temp.5gene.temp),1]
  adar.rna.cor.mat.p[,i] <- adar.rna.cor[[i]][unique(edit.temp.5gene.temp),2]
  adar2.rna.cor.mat.p[,i] <- adar2.rna.cor[[i]][unique(edit.temp.5gene.temp),2]
  adar12.cor[i] <- adar.rna.cor[[i]]['ADARB1',1]
}

#remove intergenic edits
edits.union.temp <- edits.union[nchar(edit.temp.5gene) < 20]
edits.pvalues.temp <- edits.pvalues[nchar(edit.temp.5gene) < 20,]
edits.freq.temp <- edits.freq[nchar(edit.temp.5gene) < 20,]
edits.mat.temp <- edits.mat[nchar(edit.temp.5gene) < 20,]
edit.temp.5gene.temp <- edit.temp.5gene[nchar(edit.temp.5gene) < 20]

#create edit--ADAR1/2 RNA correlation lists
adar.edit.rna.cor <- matrix(NA,nrow=length(edits.union.temp),ncol=14)
adar2.edit.rna.cor <- matrix(NA,nrow=length(edits.union.temp),ncol=14)
rownames(adar.edit.rna.cor) <- edits.union.temp
rownames(adar2.edit.rna.cor) <- edits.union.temp
for(j in 1:14){
  #remove intergenic edits
  adar.edit.rna.cor[intersect(rownames(adar.edit.cor[[j]][-grep('intergenic',rownames(adar.edit.cor[[j]])),]),edits.union.temp),j] <- 
    adar.edit.cor[[j]][-grep('intergenic',rownames(adar.edit.cor[[j]])),][intersect(rownames(adar.edit.cor[[j]][-grep('intergenic',rownames(adar.edit.cor[[j]])),]),edits.union.temp),2]
  adar2.edit.rna.cor[intersect(rownames(adar.edit.cor[[j]][-grep('intergenic',rownames(adar.edit.cor[[j]])),]),edits.union.temp),j] <- 
    adar2.edit.cor[[j]][-grep('intergenic',rownames(adar2.edit.cor[[j]])),][intersect(rownames(adar.edit.cor[[j]][-grep('intergenic',rownames(adar.edit.cor[[j]])),]),edits.union.temp),2]
}

#filter based on significance of edit-target RNA correlation
edits.pvalues.rowsums <- rowSums(edits.pvalues.temp < 0.1,na.rm=TRUE)
edits.mat.temp[edits.pvalues.temp >= 0.1] = 0
edits.rowmeans <- vector(mode='numeric',dim(edits.mat.temp)[1])
edits.rowmeans.sign <- matrix(0,nrow=dim(edits.mat.temp)[1],ncol=2)
for(i in 1:dim(edits.mat.temp)[1]){
  edits.rowmeans[i] <- sum(edits.mat.temp[i,][edits.pvalues.temp[i,] < 0.1],na.rm = TRUE)/sum(edits.pvalues.temp[i,] < 0.1,na.rm=TRUE)
  edits.rowmeans.sign[i,] <- c(sum(edits.mat.temp[i,][edits.pvalues.temp[i,] < 0.1] > 0,na.rm = TRUE),sum(edits.mat.temp[i,][edits.pvalues.temp[i,] < 0.1] < 0,na.rm = TRUE))
}
names(edits.rowmeans) <- rownames(edits.mat.temp)
rownames(edits.rowmeans.sign) <- rownames(edits.mat.temp)
edits.rowmeans[is.na(edits.rowmeans)] <- 0
edits.pvalues.rowsums[edits.rowmeans == 0] <- 0

#filter based on edits being correlated (edits - target RNA) in the same direction
edit.temp.5gene.temp <- edit.temp.5gene.temp[(edits.rowmeans.sign[,1] > 0 & edits.rowmeans.sign[,2] == 0) | (edits.rowmeans.sign[,1] == 0 & edits.rowmeans.sign[,2] > 0) | (edits.rowmeans.sign[,1] == 0 & edits.rowmeans.sign[,2] == 0)]
edits.pvalues.rowsums <- edits.pvalues.rowsums[(edits.rowmeans.sign[,1] > 0 & edits.rowmeans.sign[,2] == 0) | (edits.rowmeans.sign[,1] == 0 & edits.rowmeans.sign[,2] > 0) | (edits.rowmeans.sign[,1] == 0 & edits.rowmeans.sign[,2] == 0)]
edits.rowmeans <- edits.rowmeans[(edits.rowmeans.sign[,1] > 0 & edits.rowmeans.sign[,2] == 0) | (edits.rowmeans.sign[,1] == 0 & edits.rowmeans.sign[,2] > 0) | (edits.rowmeans.sign[,1] == 0 & edits.rowmeans.sign[,2] == 0)]

#filter based on editing type, UTR3, UTR5, synonymous, etc...
edit.temp.5gene.temp <- edit.temp.5gene.temp[grep('UTR3',names(edits.rowmeans))]
edits.pvalues.rowsums <- edits.pvalues.rowsums[grep('UTR3',names(edits.rowmeans))]
edits.rowmeans <- edits.rowmeans[grep('UTR3',names(edits.rowmeans))]

#initialize outputs
edits.sig.temp <- list()
edits.freq.temp.list <- list()
edits.adar.cor.temp <- list()
adar.edit.sig.temp <- list()
adar2.edit.sig.temp <- list()
edits.sigPerGene <- matrix(0,nrow=length(unique(edit.temp.5gene.temp)),ncol=14)
rownames(edits.sigPerGene) <- unique(edit.temp.5gene.temp)
edits.sigPerGene.cor <- edits.sigPerGene

#convert edit-level data to gene-level lists
for(i in 1:length(unique(edit.temp.5gene.temp))){
  adar.edit.sig.temp[[i]] <- adar.edit.rna.cor[which(edit.temp.5gene.temp == unique(edit.temp.5gene.temp)[i]),]
  adar2.edit.sig.temp[[i]] <- adar2.edit.rna.cor[which(edit.temp.5gene.temp == unique(edit.temp.5gene.temp)[i]),]
  #edits.sig.temp[[i]] <- cbind(edits.pvalues.rowsums,edits.rowmeans)[which(edit.temp.5gene.temp == unique(edit.temp.5gene.temp)[i]),]
  #edits.sig.temp[[i]] <- edits.pvalues.temp[which(edit.temp.5gene.temp == unique(edit.temp.5gene.temp)[i]),]
  edits.freq.temp.list[[i]] <- edits.freq.temp[which(edit.temp.5gene.temp == unique(edit.temp.5gene.temp)[i]),]
  #edits.sigPerGene[i,u] <- sum(rowSums(edits.pvalues.temp < 0.1,na.rm=TRUE)[which(edit.temp.5gene.temp == unique(edit.temp.5gene)[i])] > u-1)
  if(i %% 1000 == 0){print(i)}
}


#convert gene level lists to gene vectors and matrices
edits.freq.cancer <- matrix(0,nrow=length(unique(edit.temp.5gene.temp)),ncol=14)
rownames(edits.freq.cancer) <- unique(edit.temp.5gene.temp)
adar.edit.cancer <- matrix(0,nrow=length(unique(edit.temp.5gene.temp)),ncol=14)
rownames(adar.edit.cancer) <- unique(edit.temp.5gene.temp)
adar2.edit.cancer <- matrix(0,nrow=length(unique(edit.temp.5gene.temp)),ncol=14)
rownames(adar2.edit.cancer) <- unique(edit.temp.5gene.temp)
#for(u in 1:14){ 
  for(i in 1:length(unique(edit.temp.5gene.temp))){
    if(length(edits.freq.temp.list[[i]]) > 14){
      edits.freq.cancer[i,] <- colMaxs(edits.freq.temp.list[[i]],na.rm=TRUE)
      adar.edit.cancer[i,] <- colMins(adar.edit.sig.temp[[i]],na.rm=TRUE)
      adar2.edit.cancer[i,] <- colMins(adar2.edit.sig.temp[[i]],na.rm=TRUE)
      #edits.sigPerGene[i,] <- colSums(edits.sig.temp[[i]] < 0.1,na.rm=TRUE)
    
    } else {
      adar.edit.cancer[i,] <- adar.edit.sig.temp[[i]]
      adar2.edit.cancer[i,] <- adar2.edit.sig.temp[[i]]
      edits.freq.cancer[i,] <- edits.freq.temp.list[[i]]
      #edits.sigPerGene[i,] <- as.numeric(edits.sig.temp[[i]] < 0.1)
    }
    #if(length(edits.sig.temp[[i]]) > 2){
      #edits.sigPerGene[i,u] <- sum(edits.sig.temp[[i]][,1] > u-1)
      #edits.sigPerGene.cor[i,u] <- median(edits.sig.temp[[i]][,2][edits.sig.temp[[i]][,2] != 0 & edits.sig.temp[[i]][,1] > u-1])
    #} else if(length(edits.sig.temp[[i]]) > 1){
      #edits.sigPerGene[i,u] <- sum(edits.sig.temp[[i]][1] > u-1)
    #}
  }
#  print(u)
#}
adar.edit.cancer[is.finite(adar.edit.cancer) == FALSE] <- NA
adar2.edit.cancer[is.finite(adar.edit.cancer) == FALSE] <- NA

#finally, use your outputs to filter your gene list
ad1 <- unique(edit.temp.5gene.temp)[rowSums(adar.rna.cor.mat.p < 0.1 & edits.freq.cancer > 0.1 & adar.edit.cancer < 0.1,na.rm=TRUE) > 7]
ad2 <- unique(edit.temp.5gene.temp)[rowSums(adar2.rna.cor.mat.p < 0.1 & edits.freq.cancer > 0.1 & adar2.edit.cancer < 0.1,na.rm=TRUE) > 7]

#annotate each gene by the number of frequently edited sites from each type
#types  <- types[-c(3,10)]
ad1.types <- matrix(0,nrow=length(ad1),ncol=length(types))
rownames(ad1.types) <- ad1 
colnames(ad1.types) <- types
for(i in 1:length(ad1)){
  for(j in 1:length(types)){
    temp52 <- intersect(grep(types[j],rownames(edits.pvalues.temp)),which(edit.temp.5gene.temp == ad1[i]))
    if(length(temp52) > 1){
      ad1.types[i,j] <- sum(rowSums(edits.freq.temp[temp52,] > 0.1 & edits.pvalues.temp[temp52,] < 0.1) > 0)
    } else if(length(temp52) > 0){
      ad1.types[i,j] <- sum(sum(edits.freq.temp[temp52,] > 0.1 & edits.pvalues.temp[temp52,] < 0.1) > 0)
    }
  }
}

ad2.types <- matrix(0,nrow=length(ad2),ncol=length(types))
rownames(ad2.types) <- ad2 
colnames(ad2.types) <- types
for(i in 1:length(ad2)){
  for(j in 1:length(types)){
    temp52 <- intersect(grep(types[j],rownames(edits.pvalues.temp)),which(edit.temp.5gene.temp == ad2[i]))
    if(length(temp52) > 1){
      ad2.types[i,j] <- sum(rowSums(edits.freq.temp[temp52,] > 0.1 & edits.pvalues.temp[temp52,] < 0.1) > 0)
    } else if(length(temp52) > 0){
      ad2.types[i,j] <- sum(sum(edits.freq.temp[temp52,] > 0.1 & edits.pvalues.temp[temp52,] < 0.1) > 0)
    }
  }
}

#create bed file of segments of genes with a bunch of significant edits
#split the edits by type, then for each type with > 0 edits, find the first and last edit and add -20 & +20 to the position
#record the chromosome and strand of the edits
#save as data.frame
ad1 <- ad1[ad1.types[,2] >= 10]
type <- types[2]
k = 0
bed.save.temp <- matrix(NA,nrow=length(ad1)*15,ncol=5)
for(i in 1:length(ad1)){
  temp52 <- intersect(grep(type,rownames(edits.pvalues.temp)),which(edit.temp.5gene.temp == ad1[i]))
  u = 1
  v = 1
  #bed.save.temp.3 <- vector(mode='numeric',length(temp52))
  #bed.save.temp.4 <- vector(mode='numeric',length(temp52)-1)
  for(j in 1:length(temp52)){
    u = u+1
    if(j < (length(temp52)-1)){
      bed.save.temp.2 <- as.numeric(strsplit(rownames(edits.freq.temp)[temp52],'\\|')[[j+1]][2])-as.numeric(strsplit(rownames(edits.freq.temp)[temp52],'\\|')[[j]][2])
    }
    #bed.save.temp.3[j+1] <- bed.save.temp.2+bed.save.temp.3[j]
    #bed.save.temp.4[j] <- bed.save.temp.2
    if(bed.save.temp.2 > 100 | j == length(temp52)){
      if(u > 5){
        k = k+1
        bed.save.temp[k,3] <- as.numeric(strsplit(rownames(edits.freq.temp)[temp52],'\\|')[[v]][2])-20
        bed.save.temp[k,c(2,5)] <- strsplit(rownames(edits.freq.temp)[temp52],'\\|')[[1]][c(1,5)]
        bed.save.temp[k,1] <- strsplit(rownames(edits.freq.temp)[temp52],'\\|')[[1]][4]
        bed.save.temp[k,4] <- as.numeric(strsplit(rownames(edits.freq.temp)[temp52],'\\|')[[j]][2])+20
      }
      v = j+1
      u = 1
    }
    #print(c(j,k,u,v,bed.save.temp[k,],bed.save.temp.2))
  }
}
bed.save.temp <- bed.save.temp[1:k,]
#output the bed file
write.table(bed.save.temp[,2:5],file='/Users/michaelsharpnack/Desktop/edits.bed',row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)

#output genes correlated with ADAR1 and ADAR2
write.csv(ad1,file='/Users/michaelsharpnack/Desktop/adar1.cor.freq.csv')
write.csv(ad2,file='/Users/michaelsharpnack/Desktop/adar2.cor.freq.csv')
write.csv(intersect(ad1,ad2),file='/Users/michaelsharpnack/Desktop/adar12.cor.freq.csv')

#read edited region - RBP motif spreadsheet and turn it into a frequency matrix
fimo.edit <- fread('/Users/michaelsharpnack/Downloads/fimo_editedRegion_gene.txt')
gene.rbp.freq <- matrix(0,nrow=length(unique(fimo.edit$`sequence name`)),ncol=length(unique(fimo.edit$`tf name`)))
rownames(gene.rbp.freq) <- sort(unique(fimo.edit$`sequence name`))
colnames(gene.rbp.freq) <- unique(fimo.edit$`tf name`)
for(i in 1:dim(gene.rbp.freq)[1]){
  rbps.temp <- fimo.edit$`tf name`[which(fimo.edit$`sequence name` == rownames(gene.rbp.freq)[i])]
  for(j in 1:length(rbps.temp)){
    gene.rbp.freq[i,rbps.temp[j]] <- gene.rbp.freq[i,rbps.temp[j]]+1
  }
}

thresh=3
gene.rbp.freq.above <- matrix(NA,nrow=0,ncol=3)
for(i in 1:dim(gene.rbp.freq)[1]){
  if(length(gene.rbp.freq[i,][gene.rbp.freq[i,] >= thresh]) > 0){
    rbps.temp <- vector(mode='character',0)
    above.temp <- cbind(rep(rownames(gene.rbp.freq)[i],length(gene.rbp.freq[i,][gene.rbp.freq[i,] >= thresh])),
                        colnames(gene.rbp.freq)[gene.rbp.freq[i,] >= thresh],
                        gene.rbp.freq[i,][gene.rbp.freq[i,] >= thresh])
    gene.rbp.freq.above <- rbind(gene.rbp.freq.above , above.temp)
  }
}
rownames(gene.rbp.freq.above) <- NULL
rbps.temp <- vector(mode='character',dim(gene.rbp.freq.above)[1])
for(i in 1:dim(gene.rbp.freq.above)[1]){
  rbps.temp[i] <- bed.save.temp[,1][which(paste(bed.save.temp[,2],':',bed.save.temp[,3],'-',bed.save.temp[,4],'(',bed.save.temp[,5],')',sep="") == gene.rbp.freq.above[i,1])]
}

rbps.gene.mat <- cbind(gene.rbp.freq.above,adar.rna.cor.mat[rbps.temp,])[order(gene.rbp.freq.above[,2]),]
colnames(rbps.gene.mat) <- c('Genomic Position','RNA binding protein','RBP sites on segment',cancers)
write.csv(rbps.gene.mat,file='/Users/michaelsharpnack/Desktop/rbps.gene.mat.csv')

#################################################################################################################################################




#plot pca of cancers with significant edits
edits.pvalues.logical <- edits.pvalues < 0.1
class(edits.pvalues.logical) <- 'numeric'
edits.pvalues.logical <- edits.pvalues.logical[rowSums(edits.pvalues.logical,na.rm=TRUE) > 0,]
ir.pca <- prcomp(t(na.omit(edits.pvalues.logical)))
library(ggbiplot)
g <- ggbiplot(ir.pca, obs.scale = 0, var.scale = 0, 
              ellipse = FALSE, groups = c(rep('1',12),'2',rep('3',1)),
              circle = FALSE,var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'right')
print(g)



adar.rna.cor <- list()
adar2.rna.cor <- list()
adar.edit.cor <- list()
adar2.edit.cor <- list()
for(j in 1:length(cancers)){
  
  print(paste("Reading in editing file for ",cancers[j]))
  temp.edit <- fread(paste('tcga_edits/',files.edit[j],sep=""))
  temp.edit.rownames <- temp.edit[[1]]
  temp.edit <- as.matrix(temp.edit[,-1])
  temp.edit <- temp.edit[,-dim(temp.edit)[2]]
  class(temp.edit) <- 'numeric'
  rownames(temp.edit) <- temp.edit.rownames
  colnames(temp.edit) <- substr(colnames(temp.edit),8+nchar(cancers[j]),23)
  
  print(paste("Reading in RNAseq file for ",cancers[j]))
  temp.rna <- fread(paste('tcga.rnaseq/',files.rna[j],sep=""))
  temp.rna <- as.matrix(temp.rna[,-1])
  temp.rna <- temp.rna[-c(1:30),]
  class(temp.rna) <- 'numeric'
  rownames(temp.rna) <- rna.geneid
  colnames(temp.rna) <- substr(colnames(temp.rna),1,12)
  
  temp.rna <- temp.rna[,intersect(colnames(temp.rna),colnames(temp.edit))]
  temp.edit <- temp.edit[,intersect(colnames(temp.rna),colnames(temp.edit))]
  
  print(paste(cancers[j],dim(temp.rna)[2]))
  
  temp.adar.rna.cor <- matrix(NA,nrow=dim(temp.rna)[1],ncol=2)
  temp.adar2.rna.cor <- matrix(NA,nrow=dim(temp.rna)[1],ncol=2)
  temp.adar.edit.cor <- matrix(NA,nrow=dim(temp.edit)[1],ncol=2)
  temp.adar2.edit.cor <- matrix(NA,nrow=dim(temp.edit)[1],ncol=2)
  
  print('calculating ADAR 1/2 -- editing frequency correlation')
  for(i in 1:dim(temp.edit)[1]){
    temp1 <- try(cor.test(temp.rna['ADAR',][is.na(temp.edit[i,]) == FALSE],temp.edit[i,][is.na(temp.edit[i,]) == FALSE],method='spearman'),silent=TRUE)
    temp2 <- try(cor.test(temp.rna['ADARB1',][is.na(temp.edit[i,]) == FALSE],temp.edit[i,][is.na(temp.edit[i,]) == FALSE],method='spearman'),silent=TRUE)
    if(length(temp1) > 1){
      temp.adar.edit.cor[i,1] <- temp1$estimate
      temp.adar.edit.cor[i,2] <- temp1$p.value
    }
    if(length(temp2) > 1){
      temp.adar2.edit.cor[i,1] <- temp2$estimate
      temp.adar2.edit.cor[i,2] <- temp2$p.value
    }
    if(i %% 1000 == 0){print(i)}
  }
  rownames(temp.adar.edit.cor) <- rownames(temp.edit)
  rownames(temp.adar2.edit.cor) <- rownames(temp.edit)
  
  print('calculating ADAR 1/2 -- RNA expression correlation')
  for(i in 1:dim(temp.rna)[1]){
    temp1 <- try(cor.test(temp.rna['ADAR',],temp.rna[i,],method='spearman'),silent=TRUE)
    temp2 <- try(cor.test(temp.rna['ADARB1',],temp.rna[i,],method='spearman'),silent=TRUE)
    if(length(temp1) > 1){
      temp.adar.rna.cor[i,1] <- temp1$estimate
      temp.adar.rna.cor[i,2] <- temp1$p.value
    }
    if(length(temp2) > 1){
      temp.adar2.rna.cor[i,1] <- temp2$estimate
      temp.adar2.rna.cor[i,2] <- temp2$p.value
    }
    if(i %% 1000 == 0){print(i)}
  }
  
  adar.rna.cor[[j]] <- temp.adar.rna.cor
  adar2.rna.cor[[j]] <- temp.adar2.rna.cor
  adar.edit.cor[[j]] <- temp.adar.edit.cor
  adar2.edit.cor[[j]] <- temp.adar2.edit.cor
}

for(i in 1:14){
  #adar.edit.cor[[i]][,2] <- p.adjust(adar.edit.cor[[i]][,2],method='BH')
  #adar2.edit.cor[[i]][,2] <- p.adjust(adar2.edit.cor[[i]][,2],method='BH')
  #adar.rna.cor[[i]][,2] <- p.adjust(adar.rna.cor[[i]][,2],method='BH')
  #adar2.rna.cor[[i]][,2] <- p.adjust(adar2.rna.cor[[i]][,2],method='BH')
  
  rownames(adar.rna.cor[[i]]) <- rna.geneid
  rownames(adar2.rna.cor[[i]]) <- rna.geneid
}

#next pull out all of the correlation data in the adar.rna.cor type objects for genes in the sigPerGene matrix
adar.cor.sigPerGene <- matrix(0,nrow=dim(edits.sigPerGene)[1],ncol=14)
rownames(adar.cor.sigPerGene) <- rownames(edits.sigPerGene)
adar2.cor.sigPerGene <- matrix(0,nrow=dim(edits.sigPerGene)[1],ncol=14)
rownames(adar2.cor.sigPerGene) <- rownames(edits.sigPerGene)
for(i in 1:14){
  adar.cor.sigPerGene[,i] <- adar.rna.cor[[i]][rownames(edits.sigPerGene),1]
  adar2.cor.sigPerGene[,i] <- adar2.rna.cor[[i]][rownames(edits.sigPerGene),1]
}

#Question: Are genes with edits more highly correlated than you would expect? Are genes with significant edits more highly correlated than you would expect?
#is this true across cancers or is it cancer dependent

hist(adar.rna.cor[[1]][setdiff(rownames(adar.rna.cor[[1]]),rownames(edits.sigPerGene)),1],col='skyblue',border=F,50)
hist(adar.cor.sigPerGene[,1],add=T,col=scales::alpha('red',.5),border=F,25)

hist(adar.rna.cor[[1]][setdiff(rownames(adar.rna.cor[[1]]),rownames(edits.sigPerGene)),1],col='skyblue',border=F,50)
hist(adar.cor.sigPerGene[rowSums(edits.sigPerGene,na.rm = TRUE) > 1,1],add=T,col=scales::alpha('red',.5),border=F,25)

hist(adar2.rna.cor[[1]][setdiff(rownames(adar2.rna.cor[[1]]),rownames(edits.sigPerGene)),1],col='skyblue',border=F,50)
hist(adar2.cor.sigPerGene[rowSums(edits.sigPerGene,na.rm = TRUE) > 1,1],add=T,col=scales::alpha('red',.5),border=F,50)


hist(adar2.cor.sigPerGene[rowSums(edits.sigPerGene,na.rm = TRUE) > 7,10],ylim=c(0,100),col='skyblue',border=F,50)
hist(adar.cor.sigPerGene[rowSums(edits.sigPerGene,na.rm = TRUE) > 7,10],add=T,col=scales::alpha('red',.5),border=F,50)

plot(adar.cor.sigPerGene[rownames(edits.sigPerGene[edits.sigPerGene[,7] > 0,][rowSums(edits.sigPerGene[edits.sigPerGene[,7] > 0,]) > 15,]),2],
  adar2.cor.sigPerGene[rownames(edits.sigPerGene[edits.sigPerGene[,7] > 0,][rowSums(edits.sigPerGene[edits.sigPerGene[,7] > 0,]) > 15,]),2],
  pch=19,col='skyblue',cex=2,xlab='ADAR-cor',ylab='ADAR2-cor')
  abline(h=0,v=0,lwd=3,col='red')

plot(adar.cor.sigPerGene[rownames(edits.sigPerGene[edits.sigPerGene[,7] > 0,][rowSums(edits.sigPerGene[edits.sigPerGene[,7] > 0,]) > 15,]),10]*
     edits.sigPerGene.cor[edits.sigPerGene[,7] > 0,][rowSums(edits.sigPerGene[edits.sigPerGene[,7] > 0,]) > 15,7],
     adar2.cor.sigPerGene[rownames(edits.sigPerGene[edits.sigPerGene[,7] > 0,][rowSums(edits.sigPerGene[edits.sigPerGene[,7] > 0,]) > 15,]),10]*
     edits.sigPerGene.cor[edits.sigPerGene[,7] > 0,][rowSums(edits.sigPerGene[edits.sigPerGene[,7] > 0,]) > 15,7],
     pch=19,col='skyblue',cex=2,xlab='ADAR-cor',ylab='ADAR2-cor')
     abline(h=0,v=0,lwd=3,col='red')

################################################################################################################################################################################
#Validate editing correlation in lung adenocarcinoma TCGA with your other 2 datasets + ccle luad cell lines
#Q1: what are the most edited sites in all datasets
#Q2: what is the correlation between ADAR1 & ADAR2 in each dataset
#Q3: which edits are correlated with ADAR1 or 2 expression in each dataset
#Q4: which edits are correlated with their gene's expression
#Q5: what is the correlation between ADAR each gene
#Q6: are any edits predictive of survival?
#Q7: which edits are significantly different b/w tumor/normal in Carbone early stage + tcga?
#Q8: Annotate the genes thare are significant in Q7 and check to see if there is differential RNAseq expression of these genes as well in tumor/normal
#Q9: Add in information about copy number and tumor immune infiltrate levels
     

#load in tcga luad rnaseq + editing data
{
cancers <- c('BLCA','BRCA','CESC','CRCC','GBM','HNSC','KICH','KIRC','KIRP','LGG','LIHC','LUAD','LUSC','PRAD','STAD','THCA','UCEC')
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/tcga_edits/')
files.edit <- dir()
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/tcga.rnaseq/')
files.rna <- dir()
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query')
temp.rna.geneid <- read.csv(paste('tcga.rnaseq/',files.rna[1],sep=""),stringsAsFactors=FALSE,sep='\t')[,1]
temp.names <- strsplit(temp.rna.geneid,"\\|")
temp.names <- temp.names[-c(1:30)]
temp.names.2 <- vector(mode='numeric',length(temp.names))
for(i in 1:length(temp.names)){
  temp.names.2[i] <- temp.names[[i]][1]
}
rna.geneid <- temp.names.2
rm(temp.names,temp.names.2,temp.rna.geneid)
temp.edit <- fread(paste('tcga_edits/',files.edit[12],sep=""))
temp.edit.rownames <- temp.edit[[1]]
temp.edit <- as.matrix(temp.edit[,-1])
temp.edit <- temp.edit[,-dim(temp.edit)[2]]
class(temp.edit) <- 'numeric'
rownames(temp.edit) <- temp.edit.rownames
temp.edit.normal <- temp.edit[,grep('Normal',colnames(temp.edit))]
temp.edit <- temp.edit[,grep('Tumor',colnames(temp.edit))]
colnames(temp.edit) <- substr(colnames(temp.edit),8+nchar(cancers[12]),23)
colnames(temp.edit.normal) <- substr(colnames(temp.edit.normal),9+nchar(cancers[12]),24)

print(paste("Reading in RNAseq file for ",cancers[12]))
temp.rna <- fread(paste('tcga.rnaseq/',files.rna[12],sep=""))
temp.rna <- as.matrix(temp.rna[,-1])
temp.rna <- temp.rna[-c(1:30),]
class(temp.rna) <- 'numeric'
rownames(temp.rna) <- rna.geneid
temp.rna.normal <- temp.rna[,substr(colnames(temp.rna),14,14) == '1']
temp.rna <- temp.rna[,substr(colnames(temp.rna),14,14) == '0']
colnames(temp.rna.normal) <- substr(colnames(temp.rna.normal),1,12)
colnames(temp.rna) <- substr(colnames(temp.rna),1,12)

tcga.rna <- temp.rna[,intersect(colnames(temp.rna),colnames(temp.edit))]
tcga.edit <- temp.edit[,intersect(colnames(temp.rna),colnames(temp.edit))]
rm(temp.edit,rna.geneid,temp.rna,cancers,files.edit,files.rna,i,temp.edit.rownames)
}
#load in ccle luad rnaseq + editing data
{
dir.edits.ccle <- dir('/Volumes/Mac_Backup/output/')
dir.rnaseq.ccle <- dir('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/ccle/ccle_gene_quantification/')
key.edits.ccle <- read.csv('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/ccle/gdc_manifest.2017-04-21T17_46_11.622648.tsv.txt',stringsAsFactors=FALSE,sep='\t')
rownames(key.edits.ccle) <- key.edits.ccle$id
key.edits.ccle <- key.edits.ccle[substr(dir.edits.ccle,1,36),]
meta.rnaseq.ccle <- read.csv('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/ccle/cgHub_CCLE_RNA-seq_metadata_summary.txt',stringsAsFactors=FALSE,sep='\t')
rownames(meta.rnaseq.ccle) <- meta.rnaseq.ccle$analysis_id
meta.rnaseq.ccle <- meta.rnaseq.ccle[substr(dir.rnaseq.ccle,1,36),]
rownames(meta.rnaseq.ccle) <- meta.rnaseq.ccle$filename
names(dir.edits.ccle) <- key.edits.ccle$filename
names(dir.rnaseq.ccle) <- meta.rnaseq.ccle$filename
dir.edits.ccle <- dir.edits.ccle[intersect(names(dir.edits.ccle),names(dir.rnaseq.ccle))]
dir.rnaseq.ccle <- dir.rnaseq.ccle[intersect(names(dir.edits.ccle),names(dir.rnaseq.ccle))]
meta.int <- meta.rnaseq.ccle[names(dir.rnaseq.ccle),]
rm(meta.rnaseq.ccle,key.edits.ccle)
sample.info <- fread('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/ccle/CCLE_sample_info_file_2012-10-18.txt')
ccle.names <- vector(mode = 'character',dim(meta.int)[1])
for(i in 1:dim(meta.int)[[1]]){
  ccle.names[i] = strsplit(rownames(meta.int),'\\.')[[i]][2]  
}
for(i in 1:dim(sample.info)[1]){
  if(length(strsplit(sample.info[[2]],' ')[[i]]) > 1){
    sample.info[[2]][i] <- paste(strsplit(sample.info[[2]],' ')[[i]][1],strsplit(sample.info[[2]],' ')[[i]][2],sep='_')
  }
}
sample.info[[7]][grep('A549_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('ABC1_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('HCC2935_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('NCIH1435_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('NCIH1568_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('NCIH1793_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('NCIH1838_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('NCIH1944_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('NCIH1975_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('NCIH2030_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('NCIH23_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('NCIH522_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('NCIH838_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('PC14_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('RERFLCKJ_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
sample.info[[7]][grep('RERFLCMS_LUNG',sample.info[[1]])] <- 'adenocarcinoma'
adeno <- sample.info[sample.info[[5]] == 'lung' & (sample.info[[7]] == 'adenocarcinoma' | sample.info[[7]] == 'bronchioloalveolar_adenocarcinoma'),][[2]]
meta.int.temp <- data.frame()
count = 0
for(i in 1:length(adeno)){
  if(length(which(ccle.names == adeno[i])) > 0){
    count = count+1
    meta.int.temp <- rbind(meta.int.temp,meta.int[which(ccle.names == adeno[i]),])
  }  
}
meta.int <- meta.int.temp
dir.edits.ccle <- dir.edits.ccle[rownames(meta.int)]
dir.rnaseq.ccle <- dir.rnaseq.ccle[rownames(meta.int)]
rm(meta.int.temp,sample.info,adeno,ccle.names,count,i)
setwd('/Volumes/Mac_Backup/output/')
temp <- fread(dir.edits.ccle[1],blank.lines.skip=TRUE,skip=3)
temp2 <- paste(temp[[1]],temp[[2]],sep='|')
edits.union.freq <- matrix(0,nrow=dim(tcga.edit)[1],ncol=length(dir.edits.ccle))
edits.union.split <- strsplit(rownames(tcga.edit),'\\|')
edits.union.ccle <- vector(mode='character',length(edits.union.split))
for(i in 1:length(edits.union.split)){
  edits.union.ccle[i] <- paste(edits.union.split[[i]][1],edits.union.split[[i]][2],sep='|')
}
edits.union.ccle <- gsub("chr","",edits.union.ccle)
rm(edits.union.split)
rownames(edits.union.freq) <- edits.union.ccle
for(i in 1:length(dir.edits.ccle)){
  print(i)
  temp.edit <- fread(dir.edits.ccle[i],blank.lines.skip=TRUE,skip=3)
  temp.edit.11 <- as.numeric(temp.edit[[7]])
  temp.edit.13 <- as.numeric(temp.edit[[9]])
  temp.edit.13[temp.edit.11 < 10] = NA
  names(temp.edit.13) <- temp2
  edits.union.freq[,i] <- temp.edit.13[edits.union.ccle]
}
ccle.edit <- edits.union.freq
rm(edits.union.freq,temp,temp.edit,dir.edits.ccle,edits.union.ccle,i,temp.edit.11,temp.edit.13,temp2)
ccle.names <- vector(mode='character',dim(meta.int)[1])
for(i in 1:dim(meta.int)[1]){
  ccle.names[i] <- strsplit(meta.int$filename,'\\.')[[i]][2]
}
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/ccle/ccle_gene_quantification/')
ccle.rna <- matrix(0,nrow=57997,ncol=length(dir.rnaseq.ccle))
for(i in 1:length(dir.rnaseq.ccle)){
  ccle.rna[,i] <- fread(dir.rnaseq.ccle[i])[[3]]
  if(i %% 10 == 0){print(i)}
}
rownames(ccle.rna) <- fread(dir.rnaseq.ccle[i])[[1]]
library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
bm.key <- getBM(c('ensembl_gene_id','hgnc_symbol'),values = rownames(ccle.rna), mart=mart)
ccle.rna <- ccle.rna[intersect(bm.key$ensembl_gene_id,rownames(ccle.rna)),]
bm.key <- getBM(c('ensembl_gene_id','hgnc_symbol'),values = rownames(ccle.rna), mart=mart)
bm.key <- bm.key[duplicated(bm.key[,1]) == FALSE,]
rownames(bm.key) <- bm.key[,1]
bm.key <- bm.key[nchar(bm.key[,2]) > 0,]
bm.key <- bm.key[duplicated(bm.key[,2]) == FALSE,]
bm.key <- bm.key[intersect(rownames(bm.key),rownames(ccle.rna)),]
ccle.rna <- ccle.rna[rownames(bm.key),]
rownames(ccle.rna) <- bm.key[,2]
rm(bm.key,mart)
colnames(ccle.edit) <- ccle.names
colnames(ccle.rna) <- ccle.names
rm(ccle.names,i,meta.int)
}
#load in Carbone stage IV rnaseq + editing data
{

load("~/Desktop/David/Metastasis/stageIV.03292017.RData")
library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
bm.key <- getBM(c('entrezgene','hgnc_symbol'),values = rownames(expression), mart=mart)
expression <- expression[intersect(bm.key$entrezgene,rownames(expression)),]
bm.key <- getBM(c('entrezgene','hgnc_symbol'),values = rownames(expression), mart=mart)
bm.key <- bm.key[duplicated(bm.key[,1]) == FALSE,]
bm.key <- bm.key[-1,]
rownames(bm.key) <- bm.key[,1]
bm.key <- bm.key[nchar(bm.key[,2]) > 0,]
bm.key <- bm.key[duplicated(bm.key[,2]) == FALSE,]
bm.key <- bm.key[intersect(rownames(bm.key),rownames(expression)),]
expression <- expression[rownames(bm.key),]
rownames(expression) <- bm.key[,2]
rm(bm.key,mart)
stage4.rna <- expression
rm(expression)

setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/s4_edits_frequency/')
files <- dir()
temp <- fread(files[1],blank.lines.skip=TRUE,skip=3)
temp2 <- paste(temp[[1]],temp[[2]],sep='|')
rm(temp)
edits.union.freq <- matrix(0,nrow=dim(tcga.edit)[1],ncol=length(files))
edits.union.split <- strsplit(rownames(tcga.edit),'\\|')
edits.union.ccle <- vector(mode='character',length(edits.union.split))
for(i in 1:length(edits.union.split)){
  edits.union.ccle[i] <- paste(edits.union.split[[i]][1],edits.union.split[[i]][2],sep='|')
}
rm(edits.union.split)
rownames(edits.union.freq) <- edits.union.ccle
for(i in 1:length(files)){
  print(i)
  temp.edit <- fread(files[i],blank.lines.skip=TRUE,skip=3)
  temp.edit.11 <- as.numeric(temp.edit[[7]])
  temp.edit.13 <- as.numeric(temp.edit[[9]])
  temp.edit.13[temp.edit.11 < 10] = NA
  names(temp.edit.13) <- temp2
  edits.union.freq[,i] <- temp.edit.13[edits.union.ccle]
}
stage4.edit <- edits.union.freq
stage4.rna <- stage4.rna[,-74]
colnames(stage4.edit) <- colnames(stage4.rna) 
rm(edits.union.freq,clinical.data,temp.edit.13,temp.edit.11,temp.edit,files,i,temp2)
}

#load in Carbone early stage rnaseq + editing data
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/pg.counts/')
files <- dir()
names <- fread(files[1])[[1]]
pg.normal <- matrix(0,nrow=length(names),ncol=41)
pg.tumor <- matrix(0,nrow=length(names),ncol=41)
rownames(pg.normal) <- names
colnames(pg.normal) <- seq(from=1,to=41,by=1)
rownames(pg.tumor) <- names
colnames(pg.tumor) <- seq(from=1,to=41,by=1)
for(i in 1:length(files)){
  if(i %% 2 == 1){
    pg.normal[,(i-1)/2+1] <- 1000000*as.numeric(fread(files[i])[[2]])/sum(as.numeric(fread(files[i])[[2]]))
  } else {
    pg.tumor[,i/2] <- 1000000*as.numeric(fread(files[i])[[2]])/sum(as.numeric(fread(files[i])[[2]]))
  }
}

setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/pg.edits.freq/')
files <- dir()
names <- fread(files[1],blank.lines.skip=TRUE)[[1]]
temp <- fread(files[1],blank.lines.skip=TRUE,skip=3)
temp2 <- paste(temp[[1]],temp[[2]],sep='|')
rm(temp)
edits.union.freq <- matrix(0,nrow=dim(tcga.edit)[1],ncol=length(files)/2)
edits.union.freq.normal <- matrix(0,nrow=dim(tcga.edit)[1],ncol=length(files)/2)
edits.union.split <- strsplit(rownames(tcga.edit),'\\|')
edits.union.ccle <- vector(mode='character',length(edits.union.split))
for(i in 1:length(edits.union.split)){
  edits.union.ccle[i] <- paste(edits.union.split[[i]][1],edits.union.split[[i]][2],sep='|')
}
rm(edits.union.split)
rownames(edits.union.freq) <- edits.union.ccle
rownames(edits.union.freq.normal) <- edits.union.ccle
edits.union.ccle <- gsub("chr","",edits.union.ccle)
for(i in 1:length(files)){
  print(i)
  temp.edit <- fread(files[i],blank.lines.skip=TRUE,skip=3)
  temp.edit.11 <- as.numeric(temp.edit[[7]])
  temp.edit.13 <- as.numeric(temp.edit[[9]])
  temp.edit.13[temp.edit.11 < 10] = NA
  names(temp.edit.13) <- temp2
  if(i %% 2 == 1){
    edits.union.freq.normal[,(i-1)/2+1] <- temp.edit.13[edits.union.ccle]
  } else {
    edits.union.freq[,i/2] <- temp.edit.13[edits.union.ccle]
  }
}
pg.edit <- edits.union.freq
pg.edit.normal <- edits.union.freq.normal
rownames(pg.edit) <- rownames(tcga.edit)
rownames(pg.edit.normal) <- rownames(tcga.edit)
rm(temp.edit,edits.union.ccle,edits.union.freq,edits.union.freq.normal,files,names,temp.edit.11,temp.edit.13,temp2)
     
edits <- list('tcga' = tcga.edit,'ccle' = ccle.edit,'stage4' = stage4.edit, 'pg' = pg.edit)
rna <- list('tcga' = tcga.rna,'ccle' = ccle.rna,'stage4' = stage4.rna, 'pg' = pg.tumor)
rna[[3]] <- rna[[3]]/colSums(rna[[3]])*1000000
temp <- matrix(0,nrow=dim(rna[[3]])[1],ncol=length(unique(colnames(rna[[3]]))))
temp2 <- matrix(0,nrow=dim(edits[[3]])[1],ncol=length(unique(colnames(rna[[3]]))))
for(i in 1:length(unique(colnames(rna[[3]])))){
  temp.1 <- rna[[3]][,which(colnames(rna[[3]]) == unique(colnames(rna[[3]]))[i])]
  temp.1.2 <- edits[[3]][,which(colnames(rna[[3]]) == unique(colnames(rna[[3]]))[i])]
  if(is.null(dim(temp.1))){
    temp[,i] <- temp.1
    temp2[,i] <- temp.1.2
  } else {
    temp[,i] <- temp.1[,which.max(colSums(temp.1))]
    temp2[,i] <- temp.1.2[,which.max(colSums(temp.1))]
  }
}
colnames(temp) <- unique(colnames(rna[[3]]))
colnames(temp2) <- unique(colnames(rna[[3]]))
rownames(temp) <- rownames(rna[[3]])
rna[[3]] <- temp
edits[[3]] <- temp2
rm(temp,temp2,temp.1,temp.1.2)
n.studies = 4
adar12.cor <- vector(mode='numeric',n.studies)
for(i in 1:4){
  rna[[i]] <- rna[[i]][intersect(intersect(intersect(rownames(ccle.rna),rownames(stage4.rna)),rownames(tcga.rna)),rownames(pg.tumor)),]
}
rownames(edits[[2]]) <- rownames(edits[[1]])
rownames(edits[[3]]) <- rownames(edits[[1]])
adar1.edit.cor <- matrix(0,nrow=dim(edits[[1]])[1],ncol=n.studies)
rownames(adar1.edit.cor) <- rownames(edits[[1]])
adar2.edit.cor <- matrix(0,nrow=dim(edits[[1]])[1],ncol=n.studies)
rownames(adar1.edit.cor) <- rownames(edits[[1]])
adar1.gene.cor <- matrix(0,nrow=dim(rna[[1]])[1],ncol=n.studies)
rownames(adar1.gene.cor) <- rownames(rna[[1]])
adar2.gene.cor <- matrix(0,nrow=dim(rna[[1]])[1],ncol=n.studies)
rownames(adar2.gene.cor) <- rownames(rna[[1]])
edit.rna.all <- list()
for(u in 1:n.studies){
  #ADAR1 - ADAR2 correlation
  adar12.cor[u] <- cor(rna[[u]]['ADAR',],rna[[u]]['ADARB1',],method='spearman')
  #editing frequency to ADAR1/2 correlation
  for(j in 1:dim(edits[[1]])[1]){
    adar1.edit.cor[j,u] <- cor(edits[[u]][j,][is.na(edits[[u]][j,]) == FALSE],rna[[u]]['ADAR',][is.na(edits[[u]][j,]) == FALSE],method='spearman')
    adar2.edit.cor[j,u] <- cor(edits[[u]][j,][is.na(edits[[u]][j,]) == FALSE],rna[[u]]['ADARB1',][is.na(edits[[u]][j,]) == FALSE],method='spearman')
  }
  #ADAR to gene correlation
  for(j in 1:dim(rna[[u]])[1]){
    adar1.gene.cor[j,u] <- cor(rna[[u]][j,],rna[[u]]['ADAR',],method='spearman')
    adar2.gene.cor[j,u] <- cor(rna[[u]][j,],rna[[u]]['ADARB1',],method='spearman')
  }
  
  #editing frequency to target gene correlation
  edit.rna.cor <- vector(mode='numeric',100000)
  edit.rna.cor.labels <- matrix(NA,nrow=100000,ncol=2)
  edit.rna.cor.nsamples <- vector(mode='numeric',100000)
  edit.rna.cor.mean.freq <- vector(mode='numeric',100000)
  edit.rna.cor.pvalues <- vector(mode='numeric',100000)
  k = 0
  for(i in 1:dim(edits[[u]])[1]){
    temp = strsplit(rownames(edits[[u]])[i],"\\|")[[1]]
    if(temp[3] == 'intergenic'){
      temp = strsplit(rownames(edits[[u]])[i],"\\(")[[1]]
      if(strsplit(temp[2],"\\,")[[1]][2] == "NONE"){
        if(length(which(rownames(rna[[u]]) == strsplit(temp[1],"\\|")[[1]][4])) == 1){
          k = k+1 
          edit.rna.cor.labels[k,1] = rownames(edits[[u]])[i]
          edit.rna.cor.labels[k,2] = strsplit(temp[1],"\\|")[[1]][4]
          temp10 <- try(cor.test(rna[[u]][strsplit(temp[1],"\\|")[[1]][4],is.na(edits[[u]][i,]) == FALSE],
                                 edits[[u]][i,is.na(edits[[u]][i,]) == FALSE],method='spearman'))
          if(length(temp10) > 1){
            edit.rna.cor[k] = temp10$estimate
            edit.rna.cor.pvalues[k] = temp10$p.value
          }
          edit.rna.cor.nsamples[k] <- sum(is.na(edits[[u]][i,]) == FALSE)
          edit.rna.cor.mean.freq[k] <- mean(edits[[u]][i,][is.na(edits[[u]][i,]) == FALSE])
        }
      } else{
        if(length(which(rownames(rna[[u]]) == strsplit(temp[1],"\\|")[[1]][4])) == 1){
          k = k+1
          edit.rna.cor.labels[k,1] = rownames(edits[[u]])[i]
          edit.rna.cor.labels[k,2] = strsplit(temp[1],"\\|")[[1]][4]
          temp10 <- try(cor.test(rna[[u]][strsplit(temp[1],"\\|")[[1]][4],is.na(edits[[u]][i,]) == FALSE],
                                 edits[[u]][i,is.na(edits[[u]][i,]) == FALSE],method='spearman'))
          if(length(temp10) > 1){
            edit.rna.cor[k] = temp10$estimate
            edit.rna.cor.pvalues[k] = temp10$p.value
          }
          edit.rna.cor.nsamples[k] <- sum(is.na(edits[[u]][i,]) == FALSE)
          edit.rna.cor.mean.freq[k] <- mean(edits[[u]][i,][is.na(edits[[u]][i,]) == FALSE])
          
        }
        if(length(which(rownames(rna[[u]]) == strsplit(temp[2],"\\,")[[1]][2])) == 1){
          k = k+1
          edit.rna.cor.labels[k,2] = strsplit(temp[2],"\\,")[[1]][2]
          edit.rna.cor.labels[k,1] = rownames(edits[[u]])[i]
          temp10 <- try(cor.test(rna[[u]][strsplit(temp[2],"\\,")[[1]][2],is.na(edits[[u]][i,]) == FALSE],
                                 edits[[u]][i,is.na(edits[[u]][i,]) == FALSE],method='spearman'))
          if(length(temp10) > 1){
            edit.rna.cor[k] = temp10$estimate
            edit.rna.cor.pvalues[k] = temp10$p.value
          }
          edit.rna.cor.nsamples[k] <- sum(is.na(edits[[u]][i,]) == FALSE)
          edit.rna.cor.mean.freq[k] <- mean(edits[[u]][i,][is.na(edits[[u]][i,]) == FALSE])
          
        }
      }
      
    } else {
      if(length(which(rownames(rna[[u]]) == temp[4])) == 1){
        k = k+1
        edit.rna.cor.labels[k,1] = rownames(edits[[u]])[i]
        edit.rna.cor.labels[k,2] = temp[4]
        temp10 <- try(cor.test(rna[[u]][temp[4],is.na(edits[[u]][i,]) == FALSE],
                               edits[[u]][i,is.na(edits[[u]][i,]) == FALSE],method='spearman'))
        if(length(temp10) > 1){
          edit.rna.cor[k] = temp10$estimate
          edit.rna.cor.pvalues[k] = temp10$p.value
        }
        edit.rna.cor.nsamples[k] <- sum(is.na(edits[[u]][i,]) == FALSE)
        edit.rna.cor.mean.freq[k] <- mean(edits[[u]][i,][is.na(edits[[u]][i,]) == FALSE])
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
  edit.rna.all[[u]] <- cbind(edit.rna.cor.labels,edit.rna.cor.nsamples,edit.rna.cor.mean.freq,edit.rna.cor,edit.rna.cor.pvalues) 
}

###########################################################################################################################################################################
#Tumor normal comparison - TCGA & pg

tn.names <- intersect(intersect(intersect(colnames(tcga.rna),colnames(temp.rna.normal)),colnames(tcga.edit)),colnames(temp.edit.normal))
temp.rna <- tcga.rna[,tn.names]
temp.rna.normal <- temp.rna.normal[,tn.names]
temp.edit <- tcga.edit[,tn.names]
temp.edit.normal <- temp.edit.normal[,tn.names]
temp.rna <- temp.rna[intersect(rownames(temp.rna),rownames(pg.tumor)),]
temp.rna.normal <- temp.rna.normal[intersect(rownames(temp.rna),rownames(pg.tumor)),]
pg.tumor <- pg.tumor[intersect(rownames(temp.rna),rownames(pg.tumor)),]
pg.normal <- pg.normal[intersect(rownames(temp.rna),rownames(pg.tumor)),]

#get data on each edit's types, gene, etc
edits.annot <- matrix(NA,nrow=dim(pg.edit)[1],ncol=2)
names.temp <- strsplit(rownames(pg.edit),'\\|')
for(i in 1:dim(pg.edit)[1]){
  edits.annot[i,1] <- names.temp[[i]][4]
  edits.annot[i,2] <- names.temp[[i]][3]
}

#figure 1A - TN ADAR expression
plot(c(rep(1,57),rep(2,57)),log2(c(temp.rna['ADAR',],temp.rna.normal['ADAR',])),xlim=c(0.75,2.25),pch=19,col='darkgreen')
segments(rep(1,57),log2(temp.rna['ADAR',]),rep(2,57),log2(temp.rna.normal['ADAR',]))
plot(c(rep(1,41),rep(2,41)),c(pg.tumor['ADAR',],pg.normal['ADAR',]),xlim=c(0.75,2.25),pch=19,col='skyblue')
segments(rep(1,41),pg.tumor['ADAR',],rep(2,41),pg.normal['ADAR',])

#figure 1B - ADAR expression vs mean editing frequency
plot(log2(temp.rna['ADAR',])-log2(temp.rna.normal['ADAR',]),colMeans(temp.edit,na.rm=TRUE)-colMeans(temp.edit.normal,na.rm=TRUE),pch=19,col='darkgreen')
abline(v=0,h=0,lwd=2)
plot(log2(pg.tumor['ADAR',])-log2(pg.normal['ADAR',]),colMeans(pg.edit,na.rm=TRUE)-colMeans(pg.edit.normal,na.rm=TRUE),pch=19,col='skyblue')
abline(v=0,h=0,lwd=2)

#get ttest matrix for edits & rna for both datasets

ttester <- function(data,group){
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
  return(ttest.mat)
}

tcga.edit.ttest <- ttester(cbind(temp.edit,temp.edit.normal),c(rep(1,57),rep(2,57)))
tcga.rna.ttest <- ttester(cbind(temp.rna,temp.rna.normal),c(rep(1,57),rep(2,57)))
pg.rna.ttest <- ttester(cbind(pg.tumor,pg.normal),c(rep(1,41),rep(2,41)))
pg.edit.ttest <- ttester(cbind(pg.edit,pg.edit.normal),c(rep(1,41),rep(2,41)))

#
plot(pg.edit.ttest[,1],-log10(pg.edit.ttest[,5]))
plot(tcga.edit.ttest[,1],-log10(tcga.edit.ttest[,5]))

write.csv(rownames(tcga.rna.ttest)[(tcga.rna.ttest[,5][is.na(tcga.rna.ttest[,5]) == FALSE & is.na(pg.rna.ttest[,5]) == FALSE] & pg.rna.ttest[,5][is.na(tcga.rna.ttest[,5]) == FALSE & is.na(pg.rna.ttest[,5]) == FALSE] < 0.05)],file='/Users/michaelsharpnack/Desktop/ttest.rna.intersect.csv')

#correlation of t-statistics between the two datasets
plot(pg.edit.ttest[,1],tcga.edit.ttest[,1])
cor(pg.edit.ttest[,1][is.na(pg.edit.ttest[,1]) == FALSE & is.na(tcga.edit.ttest[,1]) == FALSE],tcga.edit.ttest[,1][is.na(pg.edit.ttest[,1]) == FALSE & is.na(tcga.edit.ttest[,1]) == FALSE])
plot(pg.rna.ttest[,1],tcga.rna.ttest[,1])
cor(pg.rna.ttest[,1][is.na(pg.rna.ttest[,1]) == FALSE & is.na(tcga.rna.ttest[,1]) == FALSE],tcga.rna.ttest[,1][is.na(pg.rna.ttest[,1]) == FALSE & is.na(tcga.rna.ttest[,1]) == FALSE])

#correlate editing frequencies with ADAR1/2 RNA abundance
adar1.edit.cor <- matrix(0,nrow=dim(pg.edit)[1],ncol=4)
adar2.edit.cor <- matrix(0,nrow=dim(pg.edit)[1],ncol=4)
rownames(adar1.edit.cor) <- rownames(pg.edit)
rownames(adar2.edit.cor) <- rownames(pg.edit)
for(i in 1:dim(pg.edit)[1]){
  if(sum(is.na(pg.edit[i,]) == FALSE) > 10){
    adar1.edit.cor[i,3] <- cor(pg.edit[i,][is.na(pg.edit[i,]) == FALSE],pg.tumor['ADAR',][is.na(pg.edit[i,]) == FALSE],method='spearman')
    adar2.edit.cor[i,3] <- cor(pg.edit[i,][is.na(pg.edit[i,]) == FALSE],pg.tumor['ADARB1',][is.na(pg.edit[i,]) == FALSE],method='spearman')
    adar1.edit.cor[i,4] <- cor.test(pg.edit[i,][is.na(pg.edit[i,]) == FALSE],pg.tumor['ADAR',][is.na(pg.edit[i,]) == FALSE],method='spearman')$p.value
    adar2.edit.cor[i,4] <- cor.test(pg.edit[i,][is.na(pg.edit[i,]) == FALSE],pg.tumor['ADARB1',][is.na(pg.edit[i,]) == FALSE],method='spearman')$p.value
  }
  if(sum(is.na(temp.edit[i,]) == FALSE) > 10){
    adar1.edit.cor[i,1] <- cor(temp.edit[i,][is.na(temp.edit[i,]) == FALSE],temp.rna['ADAR',][is.na(temp.edit[i,]) == FALSE],method='spearman')
    adar2.edit.cor[i,1] <- cor(temp.edit[i,][is.na(temp.edit[i,]) == FALSE],temp.rna['ADARB1',][is.na(temp.edit[i,]) == FALSE],method='spearman')
    adar1.edit.cor[i,2] <- cor.test(temp.edit[i,][is.na(temp.edit[i,]) == FALSE],temp.rna['ADAR',][is.na(temp.edit[i,]) == FALSE],method='spearman')$p.value
    adar2.edit.cor[i,2] <- cor.test(temp.edit[i,][is.na(temp.edit[i,]) == FALSE],temp.rna['ADARB1',][is.na(temp.edit[i,]) == FALSE],method='spearman')$p.value
  }
}

#figure out which edits are controlled by which gene
adar1.edit.cor[adar1.edit.cor < 0] <- 0
adar2.edit.cor[adar2.edit.cor < 0] <- 0
colors.temp <- vector(mode='numeric',dim(adar1.edit.cor)[1])
colors.temp <- 1
colors.temp[sqrt(adar1.edit.cor[,1]*adar1.edit.cor[,3]) >= 0.3] <- 2
colors.temp[sqrt(adar2.edit.cor[,1]*adar2.edit.cor[,3]) >= 0.3] <- 3
colors.temp[sqrt(adar1.edit.cor[,1]*adar1.edit.cor[,3]) >= 0.3 & sqrt(adar2.edit.cor[,1]*adar2.edit.cor[,3]) >= 0.3] <- 4

plot(sqrt(adar1.edit.cor[,1]*adar1.edit.cor[,3]),sqrt(adar2.edit.cor[,1]*adar2.edit.cor[,3]),col=colors.temp,pch=19,cex=0.5)

#create a gene level matrix
#columns: Gene name, diff expression (high/low adar) (12), 
#adar1/2-expression cor (8), # sig edits total & by types (22), # adar1/2-edits sig cor (4),
#
gene.matrix <- matrix(0,nrow=dim(pg.rna)[1],ncol=18)
rownames(gene.matrix) <- rownames(pg.rna)
gene.matrix[,1:2] <- tcga.rna.ttest[,c(1,5)]
gene.matrix[,3:4] <- 
gene.matrix[,5:6] <- 
gene.matrix[,7:8] <- pg.rna.ttest[,c(1,5)]
gene.matrix[,9:10] <- 
gene.matrix[,11:12] <- 
gene.matrix[,13:16] <- adar1.edit.cor
gene.matrix[,17:20] <- adar2.edit.cor


  
  
rna.high <- (colMeans(temp.edit,na.rm=TRUE)-colMeans(temp.edit.normal,na.rm=TRUE)) > 0
rna.low <- (colMeans(temp.edit,na.rm=TRUE)-colMeans(temp.edit.normal,na.rm=TRUE)) < 0

both.diff <- rna.ttest.high[rna.ttest.high[,5] < 0.1 & rna.ttest.low[,5] > 0.3 & rna.ttest.high[,10] < 0.1,][rowSums(is.na(rna.ttest.high[rna.ttest.high[,5] < 0.1 & rna.ttest.low[,5] > 0.3 & rna.ttest.high[,10] < 0.1,])) < 1,]
rna.diff <- rna.ttest.high[rna.ttest.high[,5] < 0.1 & rna.ttest.low[,5] > 0.3 & (rna.ttest.high[,10] > 0.3 | is.na(rna.ttest.high[,10])),]


edit.valid <- cbind(edit.rna.all[[1]],edit.rna.all[[2]][,c(3:6)],edit.rna.all[[3]][,c(3:6)])[((edit.rna.all[[1]][,6] > 0 & edit.rna.all[[1]][,6] < 0.2) & (edit.rna.all[[2]][,6] > 0 & edit.rna.all[[2]][,6] < 0.2) & (edit.rna.all[[3]][,6] > 0 & edit.rna.all[[3]][,6] < 0.2)),
                                                                                    ][rowSums(is.na(cbind(edit.rna.all[[1]],edit.rna.all[[2]][,c(3:6)],edit.rna.all[[3]][,c(3:6)])[((edit.rna.all[[1]][,6] > 0 & edit.rna.all[[1]][,6] < 0.2) & (edit.rna.all[[2]][,6] > 0 & edit.rna.all[[2]][,6] < 0.2) & (edit.rna.all[[3]][,6] > 0 & edit.rna.all[[3]][,6] < 0.2)),
                                                                                                                                                                                   ])) == 0,]
gene.valid <- adar1.gene.cor[abs(adar1.gene.cor[,1]) > 0.3 & abs(adar1.gene.cor[,2]) > 0.3 & abs(adar1.gene.cor[,3]) > 0.3 & rowSums(is.na(adar1.gene.cor)) == 0 & (rowSums(sign(adar1.gene.cor)) == 3 | rowSums(sign(adar1.gene.cor)) == 0),]

nonsyn.ttest <- edit.ttest.high[grep('Nonsynonymous',rownames(edit.ttest.high)),][rowSums(is.na(edit.ttest.high[grep('Nonsynonymous',rownames(edit.ttest.high)),])) == 0,]
View(nonsyn.ttest[order(nonsyn.ttest[,5]),])

plot(log2(temp.rna['ADAR',])-log2(temp.rna.normal['ADAR',]),colMeans(temp.edit,na.rm=TRUE)-colMeans(temp.edit.normal,na.rm=TRUE),pch=19,col='skyblue')
abline(v=0,h=0)

for(i in 1:len)
graph 

######################################################################################################

load(paste("~/Desktop/David/Metastasis/editing_query/tcga.meth/adar.meth.",11,".03272017",sep=''))

intersect(substr(colnames(repeat.meth.save[[1]]),1,15),paste(colnames(temp.rna),'-01',sep=''))
intersect(substr(colnames(repeat.meth.save[[1]]),1,15),paste(colnames(temp.rna),'-11',sep=''))

repeat.meth.tumor <- repeat.meth.save[[1]][,substr(colnames(repeat.meth.save[[1]]),14,14) == '0']
repeat.meth.normal <- repeat.meth.save[[1]][,substr(colnames(repeat.meth.save[[1]]),14,14) == '1']

repeat.adar.cor <- matrix(0,nrow=dim(repeat.meth.save[[1]])[1],ncol=2)
repeat.adar <- repeat.meth.save[[1]]
repeat.adar <- repeat.adar[,substr(colnames(repeat.adar),14,14) == '0']
colnames(repeat.adar) <- substr(colnames(repeat.adar),1,12)
repeat.adar <- repeat.adar[,intersect(colnames(rna[[1]]),colnames(repeat.adar))]
repeat.adar.rna <- rna[[1]]['ADAR',intersect(colnames(rna[[1]]),colnames(repeat.adar))]
for(i in 1:length(repeat.adar.cor)){
  repeat.adar.cor[i,1] <- cor(repeat.adar[i,],repeat.adar.rna,method='spearman')
}
repeat.adar <- repeat.meth.save[[1]]
repeat.adar <- repeat.adar[,substr(colnames(repeat.adar),14,14) == '0']
colnames(repeat.adar) <- substr(colnames(repeat.adar),1,16)
repeat.adar <- repeat.adar[,intersect(names(cnv.adar.all[[12]]),colnames(repeat.adar))]
repeat.adar.cnv <- cnv.adar.all[[12]][intersect(names(cnv.adar.all[[12]]),colnames(repeat.adar))]
for(i in 1:length(repeat.adar.cor)){
  repeat.adar.cor[i,2] <- cor(repeat.adar[i,],repeat.adar.cnv)
}

repeat.ttest <- matrix(nrow=dim(repeat.meth.save[[1]])[1],ncol=5)
for(i in 1:dim(repeat.meth.save[[1]])[1]){
  ttest.temp <- t.test(repeat.meth.save[[1]][i,substr(colnames(repeat.meth.save[[1]]),14,14) == '1'],repeat.meth.save[[1]][i,substr(colnames(repeat.meth.save[[1]]),14,14) == '0'])
  repeat.ttest[i,1] <- ttest.temp$statistic
  repeat.ttest[i,2] <- ttest.temp$estimate[[1]]
  repeat.ttest[i,3] <- ttest.temp$estimate[[2]]
  repeat.ttest[i,4] <- ttest.temp$p.value
}
repeat.ttest[,5] <- p.adjust(repeat.ttest[,4],method='BH')
repeat.ttest <- cbind(repeat.meth.save[[2]],repeat.ttest)
