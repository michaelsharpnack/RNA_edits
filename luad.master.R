#file: luad.master.R
#This master script generates all of the figures and results for the LUAD RNA editing paper
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/')

##########################################################################################################################################
#cell reports paper figures
#2A
plot(as.numeric(data.cor$edit.rna.cor[,3]),-log10(as.numeric(data.cor$edit.rna.cor[,4])),pch=19,col='blue',cex=0.7)
abline(h=1,lwd=2,col='red')
#2B
hist(as.numeric(edit.cor.rna[,3]),1000)
#3A
apol1 <- data.cor$edit.rna.cor[data.cor$edit.rna.cor[,2] == 'APOL1',]
apol1.pos.temp <- strsplit(apol1[,1],'\\|')
apol1.position <- vector(mode='character',74)
for(i in 1:74){
  apol1.position[i] <- apol1.pos.temp[[i]][2]
}
apol1.position <- as.numeric(apol1.position)
plot(apol1.position,as.numeric(apol1[,5]),col='red',pch=19)

plot(apol1.position,as.numeric(apol1[,3]),col='skyblue',pch=19)
abline(h=0,lwd=2)
#supplementary table 3
write.csv(edit.cor.rna,file='/Users/michaelsharpnack/Desktop/sharpnack.SupplementaryTable3')

##########################################################################################################################################

#script 1: /Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/luad.load.R
#a. loads in datasets
#b. output are normalized tumor and normal RNA expression and editing frequency matrices

#output in RData file:
load('luad.load.RData')

#Notes: requires hard drive for CCLE editing data

##########################################################################################################################################

#script 2: /Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/
#a. Does correlation between ADAR1/2, target RNA, and RNA edits within each gene.
#b. Output are i. matrix of ADAR1/2 <-> target RNA cor
#ii. RNA edit <-> gene (list entry for each gene, with ADAR1/2 <-> edit cor included)
#iii. ADAR1/2 <-> RNA edit cor matrix

#output in RData file
load("/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/luad.editrnaadar.cor.RData")

#Notes: this script runs on one dataset at a time. To compare datasets you have to use a separate script

##########################################################################################################################################

#function 3: /Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/luad.edits.sig.R
#filters function 2's output to give you a list of possible regulatory edits, output is the edit + the sequence at that place

#output in RData file
load()

##########################################################################################################################################

#function 4: /Users/michaelsharpnack/Desktop/David/Metastasis/
#Runs a tumor-normal analysis on a dataset

#a. RNA diff t-test (editing/ADAR high vs editing/ADAR low tumors)
#b. Edits diff t-test
#c. Intersects results from script 2 and a/b to propose cancer associated regulation

#output in RData file
load()

##########################################################################################################################################
#run the pipeline
i = 12
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/')
library(data.table)
source('tcga.edit.load.R')
source('luad.editrnaadar.cor.R')
source('luad.edits.sig.R')
source('luad.ttest.R')
cancers <- c('BLCA','BRCA','CESC','CRCC','GBM','HNSC','KICH','KIRC','KIRP','LGG','LIHC','LUAD','LUSC','PRAD','STAD','THCA','UCEC')

for(i in 8:length(cancers)){
  study <- cancers[i]
  data <- tcga.loader(i)
  data.cor <- rna.edit.adar.cor(data[[2]],data[[4]])
  data.editsig <- edits.sig(data.cor[[1]],0.1)
  bed.save <- bed.creator(unique(data.editsig[[3]][,2]),data.editsig[[3]],5,data.editsig[[1]])
  data.ttest <- tumornormal(data[[2]][,data[[6]]],data[[3]][,data[[6]]],data[[4]][,data[[6]]],data[[5]][,data[[6]]],length(data[[6]]),study)
  #edit.cor.rna <- cbind(data.editsig[[3]],data.ttest$tt.edit[rownames(data.editsig[[3]]),],data.ttest$tt.rna[data.editsig[[3]][,2],])
  #edit.cor.rna.filtered <- edit.cor.rna[edit.cor.rna[,19] < 0.1 & edit.cor.rna[,24] < 0.1,][rowSums(is.na(edit.cor.rna[edit.cor.rna[,19] < 0.1 & edit.cor.rna[,24] < 0.1,])) == 0,]
  save.image(file = paste('adar.regulation/',cancers[i],'.RData',sep=''))
  if(dim(bed.save)[1] > 1){
    write.csv(paste(paste(bed.save[,2],bed.save[,3],sep=':'),paste(bed.save[,4],bed.save[,5],sep=':'),sep='-'),
            file=paste('rbpmap.bed/rbpmap.',cancers[i],'.csv',sep=''))
  }
}


setwd('adar.regulation/')
files <- dir()
genes.edit.reg <- matrix(0,nrow=1000,ncol=15)
genes.edit.reg.names <- vector(mode='character',1000)
k = 0
for(i in 1:15){
  load(files[i])
  if(dim(edit.cor.rna.filtered)[1] > 0){
    for(j in 1:dim(edit.cor.rna.filtered)[1]){
      if(length(which(genes.edit.reg.names == edit.cor.rna.filtered[j,2])) == 0){
        k = k+1
        genes.edit.reg.names[k] <- edit.cor.rna.filtered[j,2]
        genes.edit.reg[k,i] <- sign(as.numeric(edit.cor.rna.filtered[j,20]))
      } else {
        genes.edit.reg[which(genes.edit.reg.names == edit.cor.rna.filtered[j,2]),i] <- sign(as.numeric(edit.cor.rna.filtered[j,20]))
      }
    }
  }
}
genes.edit.reg <- genes.edit.reg[1:k,]
genes.edit.reg.names <- genes.edit.reg.names[1:k]
rownames(genes.edit.reg) <- genes.edit.reg.names
View(genes.edit.reg[order(rowSums(genes.edit.reg != 0),decreasing=TRUE),])

##########################################################################################################################################
#make plots

#Figure 2: Tumor normal analysis of TCGA LUAD
#THIS IS REALLY SUPPLEMENTARY
#2a. separating patients based on ADAR abundance/editing frequency
plot(log2(data[[2]]['ADAR',data[[6]]])-log2(data[[3]]['ADAR',data[[6]]]),colMeans(data[[4]][,data[[6]]],na.rm=TRUE)-colMeans(data[[5]][,data[[6]]],na.rm=TRUE),
     pch=19,col='skyblue')
abline(v=0,h=0,lwd=2)
abline(h=0.01,lwd=2,col='red')
plot(log2(data[[2]]['ADARB1',data[[6]]])-log2(data[[3]]['ADARB1',data[[6]]]),colMeans(data[[4]][,data[[6]]],na.rm=TRUE)-colMeans(data[[5]][,data[[6]]],na.rm=TRUE),pch=19,col='darkgreen')
abline(v=0,h=0,lwd=2)

#Just show the changes in gene expression overlayed with the changes in editing frequency in those genes
#maybe a manhattan plot with rna expression t-stat on top and editing t-stat on bot
library(qqman)
edit.to.man <- data.ttest.tcga[[1]] #[is.na(data.ttest.tcga[[1]][,5]) == FALSE,]
#edit.to.man <- data.ttest[[1]] #[is.na(data.ttest[[1]][,5]) == FALSE,]
names.temp <- strsplit(rownames(edit.to.man),'\\|')
names.chr <- vector(mode='character',length(names.temp))
names.bp <- vector(mode='numeric',length(names.temp))
for(i in 1:length(names.temp)){
  names.chr[i] <- names.temp[[i]][1]
  names.bp[i] <- as.numeric(names.temp[[i]][2])
}
names.chr <- gsub('chr','',names.chr)
names.chr[names.chr == 'X'] <- 23
names.chr[names.chr == 'Y'] <- 24
names.chr <- as.numeric(names.chr)
edits.dat <- data.frame('CHR' = names.chr,'BP' = names.bp,'P'=as.numeric(edit.to.man[,5]))
manhattan(edits.dat)



################################################################################################################
#number of significant edits in each gene in each cancer type
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/')
cancers <- c('BLCA','BRCA','CESC','CRCC','GBM','HNSC','KICH','KIRC','KIRP','LGG','LIHC','LUAD','LUSC','PRAD','STAD','THCA','UCEC')
setwd('adar.regulation/')
files <- dir()

edits.rna.cor.sig <- function(){
  edits.rna.sig <- matrix(0,nrow=length(unique(data.cor$edit.rna.cor[,2])),ncol=6)
  for(i in 1:length(unique(data.cor$edit.rna.cor[,2]))){
    temp <- which(data.cor$edit.rna.cor[,2] == unique(data.cor$edit.rna.cor[,2])[i])
    edits.rna.sig[i,1] <- sum(as.numeric(data.cor$edit.rna.cor[temp,4]) <= 0.1)
    edits.rna.sig[i,2] <- sum(as.numeric(data.cor$edit.rna.cor[temp,4]) <= 0.1 & as.numeric(data.cor$edit.rna.cor[temp,3]) < 0)
    edits.rna.sig[i,3] <- sum(as.numeric(data.cor$edit.rna.cor[temp,4]) <= 0.1 & as.numeric(data.cor$edit.rna.cor[temp,3]) > 0)
    edits.rna.sig[i,4] <- sum(as.numeric(data.cor$edit.rna.cor[temp,4]) > 0.1)
    edits.rna.sig[i,5] <- sum(as.numeric(data.cor$edit.rna.cor[temp,4]) > 0.1 & as.numeric(data.cor$edit.rna.cor[temp,3]) < 0)
    edits.rna.sig[i,6] <- sum(as.numeric(data.cor$edit.rna.cor[temp,4]) > 0.1 & as.numeric(data.cor$edit.rna.cor[temp,3]) > 0)
  }
  rownames(edits.rna.sig) <- unique(data.cor$edit.rna.cor[,2])
  return(edits.rna.sig)
}
edits.rna.sig.all <- list()
for(j in 1:17){
  print(j)
  load(files[j])
  edits.rna.sig.all[[j]] <- edits.rna.cor.sig()
}
names.temp <- vector(mode='character',0)
for(j in 1:17){
  names.temp <- union(names.temp,rownames(edits.rna.sig.all[[j]]))
}
gene.nsig <- matrix(0,nrow=length(names.temp),ncol=17)
rownames(gene.nsig) <- names.temp
for(j in 1:17){
  gene.nsig[intersect(names.temp,rownames(edits.rna.sig.all[[j]])),j] <- edits.rna.sig.all[[j]][intersect(names.temp,rownames(edits.rna.sig.all[[j]])),1]
}
colnames(gene.nsig) <- substr(files,1,4)
gene.nsig[is.na(gene.nsig)] <- 0
pheatmap(log2(gene.nsig[order(gene.nsig[,12],decreasing=TRUE)[1:50],c(10,1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17)]+1),cluster_rows = FALSE,cluster_cols=TRUE)

pheatmap(log2(gene.nsig[order(gene.nsig[,12],decreasing=TRUE)[1:50],c(12,2,1,6,8,16,10,9,11,14,13,3,5,15,4,17,7)]+1),cluster_rows = FALSE,cluster_cols=FALSE)
pheatmap(log2(gene.nsig[order(gene.nsig[,12],decreasing=TRUE)[1:50],]+1),cluster_rows = FALSE,cluster_cols=FALSE)

pheatmap(log2(gene.nsig[order(rowMeans(gene.nsig,na.rm=TRUE),decreasing=TRUE),][1:50,]+1))


################################################################################################################

#Figure 2A - enrichment of regulatory edits in different gene locations
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/')
cancers <- c('BLCA','BRCA','CESC','CRCC','GBM','HNSC','KICH','KIRC','KIRP','LGG','LIHC','LUAD','LUSC','PRAD','STAD','THCA','UCEC')
setwd('adar.regulation/')
files <- dir()
loc.types.n.all <- list()
for(j in 1:15){
  print(j)
  load(files[j])
  edit.loc.temp <- strsplit(rownames(data.cor$edit.rna.cor),'\\|')
  edit.loc <- vector(mode='character',length(edit.loc.temp))
  for(i in 1:length(edit.loc.temp)){
    edit.loc[i] <- edit.loc.temp[[i]][3]
  }
  loc.types <- sort(unique(edit.loc))
  loc.types.n <- matrix(0,nrow=length(loc.types),ncol=4)
  for(i in 1:length(loc.types)){
    loc.types.n[i,1] <- sum(data.cor$edit.rna.cor[which(edit.loc == loc.types[i]),4][is.na(data.cor$edit.rna.cor[which(edit.loc == loc.types[i]),4]) == FALSE] < 0.05)
    loc.types.n[i,2] <- length(data.cor$edit.rna.cor[which(edit.loc == loc.types[i]),4])
  }
  loc.types.n[,3] <- loc.types.n[,1]/loc.types.n[,2]
  loc.types.n[,4] <- rep(sum(loc.types.n[,1])/sum(loc.types.n[,2]),length(loc.types))
  rownames(loc.types.n) <- loc.types
  colnames(loc.types.n) <- c('sig','all','%','total %')
  loc.types.n.all[[j]] <- loc.types.n
}

loc.types.plot <- matrix(0,nrow=10,ncol=15)
colnames(loc.types.plot) <- substr(files,1,4)
rownames(loc.types.plot) <- rownames(loc.types.n.all[[1]])
for(i in 1:length(loc.types.n.all)){
  loc.types.plot[rownames(loc.types.n.all[[i]]),i] <- as.numeric(loc.types.n.all[[i]][,3])
}


#make tumor normal plots
gene <- 'APOL1'
edit <- 'chr22|36662382|UTR3|APOL1|+|Alu|no-conserve'
nona <- (is.na(data[[4]][edit,data[[6]]]) == FALSE & is.na(data[[5]][edit,data[[6]]]) == FALSE)
plot(log2(data[[2]][gene,data[[6]]][nona])-log2(data[[3]][gene,data[[6]]][nona]),
     data[[4]][edit,data[[6]]][nona]-data[[5]][edit,data[[6]]][nona],pch=19,col='blue')
abline(h=0,v=0,lwd=2)
cor(log2(data[[2]][gene,data[[6]]][nona])-log2(data[[3]][gene,data[[6]]][nona]),
     data[[4]][edit,data[[6]]][nona]-data[[5]][edit,data[[6]]][nona],method='spearman')

#make edit-ADAR/ADARB1 correlation plots
plot(log2(data[[2]]['ADAR',]),data[[4]][edit,],pch=19,col='darkgreen')
data.cor$edit.rna.cor[edit,7]
plot(log2(data[[2]]['ADARB1',]),data[[4]][edit,],pch=19,col='skyblue')
data.cor$edit.rna.cor[edit,9]


plot(log2(data[[2]]['ADAR',]),log2(data[[2]][gene,]),pch=19,col='darkgreen')
data.cor$edit.rna.cor[edit,11]
plot(log2(data[[2]]['ADARB1',]),log2(data[[2]][gene,]),pch=19,col='skyblue')
data.cor$edit.rna.cor[edit,13]

boxplot(data[[2]][gene,data[[6]]],data[[3]][gene,data[[6]]])




