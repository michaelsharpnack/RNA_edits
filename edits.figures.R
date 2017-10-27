#this script creates all of the figures in the LUAD RNA editing paper

load("~/Desktop/David/Metastasis/editing_query/adar.regulation/LUAD.RData")
library(pheatmap)
#figure 2
#A
plot(as.numeric(data.cor$edit.rna.cor[,3]),-log10(as.numeric(data.cor$edit.rna.cor[,4])),pch=19,col='blue')
abline(h=1,lwd=2,col='red')
#B
hist(as.numeric(edit.cor.rna[,3]),1000)
#C
pheatmap(log2(gene.nsig[order(gene.nsig[,12],decreasing=TRUE),][1:50,]+1),cluster_rows = FALSE,cluster_cols = FALSE)
#D
hist(cor.normal.tcga[,1],col='darkgreen',25)
abline(v=c(-.3,.3),col='purple',lwd=2)
#E
plot(cor.normal.tcga[,1],cor.normal.tcga[,4],pch=19,col='darkgreen')
abline(v=c(-.3,0,.3),lwd=2,col='purple')

#figure 3
#A
edit.cor.rna.names <- strsplit(rownames(data.cor$edit.rna.cor),'\\|')
edit.cor.position <- vector(mode='numeric',length(edit.cor.rna.names))
for(i in 1:length(edit.cor.rna.names)){
  edit.cor.position[i] <- as.numeric(edit.cor.rna.names[[i]][2])
}
plot(edit.cor.position[data.cor$edit.rna.cor[,2] == 'APOL1'],as.numeric(data.cor$edit.rna.cor[data.cor$edit.rna.cor[,2] == 'APOL1',5]),pch=19,col='red')
#B
plot(edit.cor.position[data.cor$edit.rna.cor[,2] == 'APOL1'],as.numeric(data.cor$edit.rna.cor[data.cor$edit.rna.cor[,2] == 'APOL1',3]),pch=19,col='skyblue')
#C
plot(log2(data[[2]]['ADAR',]),data[[4]][grep('36662397',rownames(data[[4]])),],pch=19,col='blue')
abline(lm(data[[4]][grep('36662397',rownames(data[[4]])),]~log2(data[[2]]['ADAR',])),lwd=2,col='red')
#D
plot(log2(data[[2]]['ADAR',]),log2(data[[2]]['APOL1',]),pch=19,col='blue')
abline(lm(log2(data[[2]]['APOL1',])~log2(data[[2]]['ADAR',])),lwd=2,col='red')
#E
plot.tempv <- data[[4]][grep('36662382',rownames(data[[4]])),data[[6]]]-data[[5]][grep('36662382',rownames(data[[4]])),data[[6]]]
plot.temph <- log2(data[[2]]['APOL1',data[[6]]]/data[[3]]['APOL1',data[[6]]])
plot(
  plot.temph,
  plot.tempv,
  pch=19,col='darkgreen',xlim=c(-2,3)
)
abline(v=0,h=0,lwd=2)
abline(lm(plot.tempv~plot.temph),lwd=2,col='purple')



#Supplementary Figure 1
ttest.edit.coding <- data.ttest$tt.edit[union(grep('Nonsyn',rownames(data.ttest$tt.edit)),grep('Syn',data.ttest$tt.edit)),]
ttest.edit.coding <- ttest.edit.coding[is.na(ttest.edit.coding[,1]) == FALSE,]
ttest.edit.noncoding <- data.ttest$tt.edit[setdiff(seq(from=1,to=dim(data.ttest$tt.edit)[1],by=1),union(grep('Nonsyn',rownames(data.ttest$tt.edit)),grep('Syn',data.ttest$tt.edit))),]
ttest.edit.noncoding <- ttest.edit.noncoding[is.na(ttest.edit.noncoding[,1]) == FALSE,]
#A
plot(as.numeric(ttest.edit.coding[,1]),-log10(as.numeric(ttest.edit.coding[,5])),pch=19)
abline(h=1,lwd=2,col='red')
#B
plot(as.numeric(ttest.edit.noncoding[,1]),-log10(as.numeric(ttest.edit.noncoding[,5])),pch=19)
abline(h=1,lwd=2,col='red')

#Supplementary Figure 2
#A
plot(as.numeric(data.cor$edit.rna.cor[,7]),-log10(as.numeric(data.cor$edit.rna.cor[,8])),pch=19)
abline(h=1,lwd=2,col='red')
#B
plot(as.numeric(data.cor$edit.rna.cor[,9]),-log10(as.numeric(data.cor$edit.rna.cor[,10])),pch=19)
abline(h=1,lwd=2,col='red')

#Supplementary Figure 3
#A
boxplot(log2(data[[2]]['ADAR',data[[6]]]),log2(data[[3]]['ADAR',data[[6]]]))
#B
boxplot(log2(data[[2]]['ADARB1',data[[6]]]),log2(data[[3]]['ADARB1',data[[6]]]))

#Supplementary Figure 4
rbPal <- colorRampPalette(c('red','blue'))
color.temp <- rbPal(10)[as.numeric(cut(as.numeric(edit.cor.rna[,3]),breaks = 10))]
edit.cor.rna.sig <- edit.cor.rna[as.numeric(edit.cor.rna[,19]) < 0.1 & as.numeric(edit.cor.rna[,24]) < 0.1,]
plot(as.numeric(edit.cor.rna.sig[,15]),as.numeric(edit.cor.rna.sig[,20]),pch=19,col=color.temp,xlim=c(-20,20),ylim=c(-20,20))
color.break.max <- vector(mode='numeric',10)
for(i in 1:10){
  color.break.max[i] <- max(as.numeric(edit.cor.rna[,3])[as.numeric(cut(as.numeric(edit.cor.rna[,3]),breaks = 10)) == i])
}
legend("topleft",title="Decile",legend=color.break.max,col =rbPal(10),pch=20,cex=0.7)

#Supplementary Figure 5
rbp.apol1 <- rbp.prediction[grep('chr22:3666',rbp.prediction[,1]),]
rbp.temp <- strsplit(substr(rbp.apol1[,2],10,nchar(rbp.apol1[,2])),"\\(")
rbp.temp1 <- vector(mode='character',length(rbp.temp))
for(i in 1:length(rbp.temp)){
  rbp.temp1[i] <- rbp.temp[[i]][1]
}
rbp.hist <- vector(mode='numeric',length(unique(rbp.temp1)))
names(rbp.hist) <- unique(rbp.temp1)
for(i in 1:length(rbp.hist)){
  rbp.hist[i] <- sum(rbp.temp1 == unique(rbp.temp1)[i])
}
barplot(sort(rbp.hist),horiz=TRUE,las=2,col='blue',cex.axis=1.5)

#Supplementary Figure 6
#A
boxplot(log2(data[[2]]['PTBP1',data[[6]]]),log2(data[[3]]['PTBP1',data[[6]]]))
#B
plot(log2(data[[2]]['ADAR',]),log2(data[[2]]['PTBP1',]),pch=19,col='darkgreen')

#Supplementary Figure 7
#A
plot(as.numeric(data.cor$edit.rna.cor[data.cor$edit.rna.cor[,2] == 'MAVS',5]),as.numeric(data.cor$edit.rna.cor[data.cor$edit.rna.cor[,2] == 'MAVS',3]),col='blue',pch=19)
#B
plot(as.numeric(data.cor$edit.rna.cor[data.cor$edit.rna.cor[,2] == 'LIMD1',5]),as.numeric(data.cor$edit.rna.cor[data.cor$edit.rna.cor[,2] == 'LIMD1',3]),col='blue',pch=19)
#C
plot(as.numeric(data.cor$edit.rna.cor[data.cor$edit.rna.cor[,2] == 'VHL',5]),as.numeric(data.cor$edit.rna.cor[data.cor$edit.rna.cor[,2] == 'VHL',3]),col='blue',pch=19)



