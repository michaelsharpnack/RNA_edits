#gives you the coordinates to edit in the RNA secondary structure fasta file
load("~/Desktop/David/Metastasis/editing_query/adar.regulation/LUAD.RData")
library(data.table)
repeatmasker <- fread('/Users/michaelsharpnack/Downloads/hg19.fa.out')
View(repeatmasker[(repeatmasker[[5]] == 'chr22') & repeatmasker[[6]] > 36662079 & repeatmasker[[7]] < 36663577,])
data.apol1 <- data.cor$edit.rna.cor[which(data.cor$edit.rna.cor[,2] == 'APOL1'),][order(data.cor$edit.rna.cor[which(data.cor$edit.rna.cor[,2] == 'APOL1'),1]),]

data.apol1.temp <- data.apol1[,1]
data.apol1.temp1 <- strsplit(data.apol1.temp,'\\|')
for(i in 1:length(data.apol1.temp1)){
  data.apol1.temp[i] <- data.apol1.temp1[[i]][2]
}
class(data.apol1.temp) <- 'numeric'

data.apol1.temp[1:31] <- data.apol1.temp[1:31]-36662141
data.apol1.temp[32:74] <- data.apol1.temp[32:74]-36662668

cor.stage4.tcga <- matrix(0,nrow=50,ncol=3) 
rownames(cor.stage4.tcga) <- rownames(gene.nsig)[order(gene.nsig[,12],decreasing=TRUE)][1:50]
for(i in 1:50){
  gene <- rownames(gene.nsig)[order(gene.nsig[,12],decreasing=TRUE)][i]
  data.gene <- data.cor$edit.rna.cor[which(data.cor$edit.rna.cor[,2] == gene),][order(data.cor$edit.rna.cor[which(data.cor$edit.rna.cor[,2] == gene),1]),]
  data.gene.stage4 <- data.cor.stage4$edit.rna.cor[which(data.cor.stage4$edit.rna.cor[,2] == gene),][order(data.cor.stage4$edit.rna.cor[which(data.cor.stage4$edit.rna.cor[,2] == gene),1]),]
  cor.stage4.tcga[i,2] <- cor.test(as.numeric(data.gene[,3]),as.numeric(data.gene.stage4[,3]))$p.value
  cor.stage4.tcga[i,1] <- cor.test(as.numeric(data.gene[,3]),as.numeric(data.gene.stage4[,3]))$estimate
}
cor.stage4.tcga[,3] <- p.adjust(cor.stage4.tcga[,2],method='BH')

cor.normal.tcga <- matrix(0,nrow=344,ncol=4) 
rownames(cor.normal.tcga) <- rownames(gene.nsig)[order(gene.nsig[,12],decreasing=TRUE)][1:344]
for(i in 1:344){
  gene <- rownames(gene.nsig)[order(gene.nsig[,12],decreasing=TRUE)][i]
  data.gene <- data.cor$edit.rna.cor[which(data.cor$edit.rna.cor[,2] == gene),][order(data.cor$edit.rna.cor[which(data.cor$edit.rna.cor[,2] == gene),1]),]
  data.gene.normal <- data.cor.normal$edit.rna.cor[which(data.cor.normal$edit.rna.cor[,2] == gene),][order(data.cor.normal$edit.rna.cor[which(data.cor.normal$edit.rna.cor[,2] == gene),1]),]
  cor.normal.tcga[i,2] <- try(cor.test(as.numeric(data.gene[,3]),as.numeric(data.gene.normal[,3]))$p.value)
  cor.normal.tcga[i,1] <- try(cor.test(as.numeric(data.gene[,3]),as.numeric(data.gene.normal[,3]))$estimate)
}
cor.normal.tcga[,3] <- p.adjust(cor.normal.tcga[,2],method='BH')
cor.normal.tcga[,4] <- sort(gene.nsig[,12],decreasing=TRUE)[1:344]


View(repeatmasker[(repeatmasker[[5]] == 'chr20') & repeatmasker[[6]] > 3850000 & repeatmasker[[7]] < 3855000,])
data.mavs <- data.mavs[-1,]

data.mavs.temp <- data.mavs[,1]
data.mavs.temp1 <- strsplit(data.mavs.temp,'\\|')
for(i in 1:length(data.mavs.temp1)){
  data.mavs.temp[i] <- data.mavs.temp1[[i]][2]
}
class(data.mavs.temp) <- 'numeric'
