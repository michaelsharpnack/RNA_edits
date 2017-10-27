#this code loads in the RNA binding protein data for signific edit regulatory segments
library(data.table)
setwd('/Users/michaelsharpnack/Desktop/David/Metastasis/editing_query/rbpmap.results/')
rbp.prediction <- matrix(NA,nrow=20000,ncol=8)
k = 0
i = 0
nlines <- system('wc -l < All_Predictions.luad.txt')
while(i <= 41337){
  line.temp <- read.csv('All_Predictions.luad.txt',skip=i,nrows=1,header=FALSE,blank.lines.skip=TRUE,sep='\t',stringsAsFactors=FALSE)
  if(length(line.temp) == 1 & substr(line.temp[[1]],1,1) == 'c'){
    segment.temp <- line.temp[[1]]
  } else if(length(line.temp) == 1 & substr(line.temp[[1]],1,1) == 'P'){
    rbp.temp <- line.temp[[1]]
  } else if(dim(line.temp)[2] == 6 && substr(line.temp[[2]],1,1) == 'c'){
    k = k+1
    rbp.prediction[k,] <- c(segment.temp,rbp.temp,as.character(line.temp))
  }
  i = i+1
  if(i %% 1000 == 0){print(i)}
}
rbp.prediction <- rbp.prediction[1:k,]

bed.loc <- paste(paste(bed.save[,2],bed.save[,3],sep=':'),paste(bed.save[,4],bed.save[,5],sep=':'),sep='-')

rbp.gene.loc <- matrix(NA,nrow=dim(rbp.prediction)[1],ncol=2)
for(i in 1:dim(rbp.prediction)[1]){
  if(length(which(bed.loc == rbp.prediction[i,1])) > 0){
    rbp.gene.loc[i,1] <- bed.save[which(bed.loc == rbp.prediction[i,1]),1]
    rbp.gene.loc[i,2] <- bed.save[which(bed.loc == rbp.prediction[i,1]),6]
  }
}
rbp.prediction <- cbind(rbp.gene.loc,rbp.prediction)

rbp.predictedgenes <- list()
rbp.predictedgenes.n <- vector(mode='numeric',0)
for(i in 1:length(unique(rbp.prediction[,4]))){
  rbp.predictedgenes[[i]] <- rbp.prediction[which(rbp.prediction[,4] == unique(rbp.prediction[,4])[i]),]
  rbp.predictedgenes.n[i] <- length(unique(rbp.prediction[which(rbp.prediction[,4] == unique(rbp.prediction[,4])[i]),1]))
}
names(rbp.predictedgenes.n) <- unique(rbp.prediction[,4])

#find number of sig edits within and without the RBP motif
edit.loc.temp <- strsplit(rownames(data.editsig[[3]]),'\\|')
edit.loc.temp.1 <- matrix(0,nrow=length(edit.loc.temp),ncol=2)
edit.loc.temp.2 <- strsplit(rbp.prediction[,3],'\\:')
edit.loc.temp.3 <- matrix(0,nrow=length(edit.loc.temp.2),3)
for(i in 1:length(edit.loc.temp)){
  edit.loc.temp.1[i,1] <- substr(edit.loc.temp[[i]][1],4,nchar(edit.loc.temp[[i]][1]))
  edit.loc.temp.1[i,2] <- edit.loc.temp[[i]][2]
}
for(i in 1:length(edit.loc.temp.2)){
  edit.loc.temp.3[i,1] <- substr(edit.loc.temp.2[[i]][1],4,nchar(edit.loc.temp.2[[i]][1]))
  edit.loc.temp.3[i,2] <- strsplit(edit.loc.temp.2[[i]][2],'-')[[1]][1]
  edit.loc.temp.3[i,3] <- strsplit(edit.loc.temp.2[[i]][2],'-')[[1]][2]
}
rbp.prediction <- cbind(rbp.prediction,matrix(0,nrow=dim(rbp.prediction)[1],ncol=4))
edit.loc.temp.4 <- strsplit(rbp.prediction[,6],'\\:')
for(i in 1:dim(rbp.prediction)[1]){
  motif.start <- as.numeric(edit.loc.temp.4[[i]][2])
  motif.end <- motif.start+nchar(rbp.prediction[i,7])-6
  edits.within.motifs <- which(as.numeric(edit.loc.temp.1[,2]) >= motif.start & as.numeric(edit.loc.temp.1[,2]) <= motif.end & edit.loc.temp.1[,1] == edit.loc.temp.3[i,1])
  edits.without.motifs <- which(as.numeric(edit.loc.temp.1[,2]) >= as.numeric(edit.loc.temp.3[i,2]) & as.numeric(edit.loc.temp.1[,2]) <= as.numeric(edit.loc.temp.3[i,3]) & edit.loc.temp.1[,1] == edit.loc.temp.3[i,1])
  rbp.prediction[i,11] <- length(edits.within.motifs)
  rbp.prediction[i,12] <- length(edits.without.motifs)
  rbp.prediction[i,13] <- sum(as.numeric(data.editsig[[3]][,3][edits.without.motifs]) > 0)
  rbp.prediction[i,14] <- sum(as.numeric(data.editsig[[3]][,3][edits.without.motifs]) < 0)
  if(i %% 100 == 0){print(i)}
}

edit.name.temp <- strsplit(rownames(data.cor$edit.rna.cor),'\\|')
edit.loc.temp.5 <- vector(mode='numeric',length(edit.name.temp))
for(i in 1:length(edit.name.temp)){
  edit.loc.temp.5[i] <- as.numeric(edit.name.temp[[i]][2])
}
edit.loc.temp.6 <- vector(mode='numeric',length(edit.loc.temp.4))
for(i in 1:length(edit.loc.temp.4)){
  edit.loc.temp.6[i] <- as.numeric(edit.loc.temp.4[[i]][2])
}

gene <- 'MAVS'
rbp <- ''
plot(edit.loc.temp.5[data.cor$edit.rna.cor[,2] == gene],data.cor$edit.rna.cor[,3][data.cor$edit.rna.cor[,2] == gene],
     pch=19,col='skyblue')
abline(h=0,lwd=2)
hist(edit.loc.temp.6[rbp.prediction[,1] == gene][grep(rbp,rbp.prediction[,4][rbp.prediction[,1] == gene])],50,xlim = c(min(edit.loc.temp.5[data.cor$edit.rna.cor[,2] == gene]),max(edit.loc.temp.5[data.cor$edit.rna.cor[,2] == gene])),col='red')
plot(edit.loc.temp.5[data.cor$edit.rna.cor[,2] == gene],data.cor$edit.rna.cor[data.cor$edit.rna.cor[,2] == gene,5],col='red',pch=19)

rbp.prediction.gene <- rbp.prediction[rbp.prediction[,1] == gene,]
rbp.prediction.gene.rbp <- vector(mode='numeric',length(unique(rbp.prediction.gene[,4])))
for(i in 1:length(unique(rbp.prediction.gene[,4]))){
  rbp.prediction.gene.rbp[i] <- sum(rbp.prediction.gene[,4] == unique(rbp.prediction.gene[,4])[i])
}
names(rbp.prediction.gene.rbp) <- substr(unique(rbp.prediction.gene[,4]),10,nchar(unique(rbp.prediction.gene[,4]))-7)
barplot(sort(rbp.prediction.gene.rbp,decreasing=TRUE)[1:25],horiz=TRUE,las=2,col='red')

range <- c(208619200,208620200)
plot(edit.loc.temp.5[data.cor$edit.rna.cor[,2] == gene],data.cor$edit.rna.cor[,3][data.cor$edit.rna.cor[,2] == gene],
     pch=19,col='skyblue',xlim=range)
abline(h=0,lwd=2)
hist(edit.loc.temp.6[rbp.prediction[,1] == gene][grep(rbp,rbp.prediction[,4][rbp.prediction[,1] == gene])],100,xlim = range,col='red')
plot(edit.loc.temp.5[data.cor$edit.rna.cor[,2] == gene],data.cor$edit.rna.cor[data.cor$edit.rna.cor[,2] == gene,5],col='red',pch=19,xlim=range)

#get distribution of RNA binding protein motifs within each gene
gene <- 'APOL1'
rbp.pergene <- vector(mode='numeric',length(unique(rbp.prediction[which(rbp.prediction[,1] == gene),4])))
names.temp <- substr(unique(rbp.prediction[which(rbp.prediction[,1] == gene),4]),9,nchar(unique(rbp.prediction[which(rbp.prediction[,1] == gene),4])))
names.temp <- strsplit(names.temp,'\\(')
for(i in 1:length(unique(rbp.prediction[which(rbp.prediction[,1] == gene),4]))){
  names(rbp.pergene)[i] <- names.temp[[i]][1]
}
for(i in 1:length(unique(rbp.prediction[which(rbp.prediction[,1] == gene),4]))){
  rbp.pergene[i] <- sum(rbp.prediction[which(rbp.prediction[,1] == gene),4] == unique(rbp.prediction[which(rbp.prediction[,1] == gene),4])[i])
}
barplot(rbp.pergene,horiz=TRUE,las=2,cex.names=0.8)











