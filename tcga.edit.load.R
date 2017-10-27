

tcga.loader <- function(cancer){
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
  
  temp.edit <- fread(paste('tcga_edits/',files.edit[cancer],sep=""))
  temp.edit.rownames <- temp.edit[[1]]
  temp.edit <- as.matrix(temp.edit[,-1])
  temp.edit <- temp.edit[,-dim(temp.edit)[2]]
  class(temp.edit) <- 'numeric'
  rownames(temp.edit) <- temp.edit.rownames
  temp.edit.normal <- temp.edit[,grep('Normal',colnames(temp.edit))]
  temp.edit <- temp.edit[,grep('Tumor',colnames(temp.edit))]
  if(nchar(cancers[cancer]) == 4){
    colnames(temp.edit) <- substr(colnames(temp.edit),8+nchar(cancers[12]),23)
    colnames(temp.edit.normal) <- substr(colnames(temp.edit.normal),9+nchar(cancers[12]),24)
  } else {
    colnames(temp.edit) <- substr(colnames(temp.edit),7+nchar(cancers[12]),23)
    colnames(temp.edit.normal) <- substr(colnames(temp.edit.normal),8+nchar(cancers[12]),24)
  }

  print(paste("Reading in RNAseq file for ",cancers[cancer]))
  temp.rna <- fread(paste('tcga.rnaseq/',files.rna[cancer],sep=""))
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
  
  tcga.names <- intersect(intersect(intersect(colnames(tcga.rna),colnames(temp.rna.normal)),colnames(tcga.edit)),colnames(temp.edit.normal))
  
  return(list(cancer,tcga.rna,temp.rna.normal,tcga.edit,temp.edit.normal,tcga.names))
}

