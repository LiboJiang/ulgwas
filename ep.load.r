ep.dat.load <- function(ep1,ep2){
  
  
  dat <- list(
    n.sum = NULL,
    sample.times = NULL,
    sex = NULL,
    pheno_table = NULL,
    geno_table = NULL,
    ind = NULL,
    exp = NULL
    )
  
  exp1f <- read.csv(ep1)
  exp2f <- read.csv(ep2)
  nc <- ncol(exp1f)
  pheno <- rbind(exp1f[4:13],exp2f[4:13])
  geno <- as.matrix(rbind(exp1f[14:nc],exp2f[14:nc]))
  
  sex <- c(as.character(exp1f[,1]),as.character(exp2f[,1]))
  
  ind <- c(exp1f[,2],exp2f[,2])
  exp <- c(exp1f[,3],exp2f[,3])
  
  geno[which(geno==9)] <- NA
  newa <- paste("a",colnames(geno),sep="")
  newd <- paste("d",colnames(geno),sep="")
  
  geno.A <- geno
  geno.D <- 1 - abs(geno.A)
  colnames(geno.A) <- newa
  colnames(geno.D) <- newd
  geno.AD <- cbind(geno.A,geno.D)
  
  nm <- ncol(geno.AD)
  AD.name <- colnames(geno.AD)
  Egeno.AD <- geno.AD
  Egname <- AD.name
  for( i in 1:(nm-1)){
    for(j in (i+1):nm){
      gtmp <- geno.AD[,i]*geno.AD[,j]
      gname <- paste(AD.name[i],AD.name[j],sep="-")
      Egeno.AD <- cbind(Egeno.AD,gtmp)
      Egname <- c(Egname,gname)
    }
  }
  colnames(Egeno.AD) <- Egname
  nan.index <- which(is.nan(colMeans(Egeno.AD,na.rm=T)))
  if(length(nan.index)>0){
    Egeno.ADN <- Egeno.AD[,-na.index]
  }else{
    Egeno.ADN <- Egeno.AD
  }
    
  dat$sample.times <- ncol(pheno)
  dat$n.sum <- nrow(pheno)
  dat$pheno_table <- pheno
  dat$geno_table <- Egeno.ADN
  dat$sex <- sex
  dat$ind <- ind
  dat$exp <- exp
  return(dat)
}

ep.sex.dat.load <- function(ep1,ep2){
  
  
  dat <- list(
    n.sum = NULL,
    sample.times = NULL,
    sex = NULL,
    pheno_table = NULL,
    geno_table = NULL,
    ind = NULL,
    exp = NULL
  )
  
  exp1f <- read.csv(ep1)
  exp2f <- read.csv(ep2)
  nc <- ncol(exp1f)
  pheno <- rbind(exp1f[4:13],exp2f[4:13])
  geno <- as.matrix(rbind(exp1f[14:nc],exp2f[14:nc]))
  
  sex <- c(as.character(exp1f[,1]),as.character(exp2f[,1]))
  ind <- c(exp1f[,2],exp2f[,2])
  exp <- c(exp1f[,3],exp2f[,3])
  
  geno[which(geno==9)] <- NA
  newa <- paste("a",colnames(geno),sep="")
  newd <- paste("d",colnames(geno),sep="")
  
  geno.A <- geno
  geno.D <- 1 - abs(geno.A)
  colnames(geno.A) <- newa
  colnames(geno.D) <- newd
  geno.AD <- cbind(geno.A,geno.D)
  
  nm <- ncol(geno.AD)
  AD.name <- colnames(geno.AD)
  Egeno.AD <- geno.AD
  Egname <- AD.name
  for( i in 1:(nm-1)){
    for(j in (i+1):nm){
      gtmp <- geno.AD[,i]*geno.AD[,j]
      gname <- paste(AD.name[i],AD.name[j],sep="-")
      Egeno.AD <- cbind(Egeno.AD,gtmp)
      Egname <- c(Egname,gname)
    }
  }
  colnames(Egeno.AD) <- Egname
  nan.index <- which(is.nan(colMeans(Egeno.AD,na.rm=T)))
  if(length(nan.index)>0){
    Egeno.ADN <- Egeno.AD[,-na.index]
  }else{
    Egeno.ADN <- Egeno.AD
  }
  sex[which(sex=="F")] <- 0 
  sex[which(sex=="M")] <- 1 
  sex <- as.integer(sex)
  s.Egeno.ADN <- sex*Egeno.ADN
  s.name <- paste(colnames(Egeno.ADN),"-s",sep="")
  sc <- c("s",Egname,s.name)
  gcno <- cbind(sex,Egeno.ADN,s.Egeno.ADN)
  colnames(gcno) <- sc
  
  
  dat$sample.times <- ncol(pheno)
  dat$n.sum <- nrow(pheno)
  dat$pheno_table <- pheno
  dat$geno_table <- gcno
  dat$sex <- sex
  dat$ind <- ind
  dat$exp <- exp
  return(dat)
}

Impute<-function(Z, impute.method){
  
  p<-dim(Z)[2]
  
  if(impute.method =="random"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-2 * maf1
      }
    }
  } else if(impute.method =="bestguess") {
    
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-round(2 * maf1)
      }
    }
    
  } else {
    stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
  }
  
  return(Z)
}

ep.imputed <- function(dat,Z,impute.method ="random"){
  
  geno <- dat$geno_table
  SNPna <- colnames(Z)
  nsnp <- dim(Z)[2]
  n.c <- unique(unlist(strsplit(SNPna,"-")))
  n.c1 <- gsub("a","",n.c)
  n.c2 <- unique(gsub("d","",n.c1))
  
  all.n <- gsub("a","",colnames(geno)[1:48])
  all.n <- gsub("d","",all.n)
  missing.index <- c()
  for(i in 1:length(n.c2)){
    dat.index <- which(n.c2[i]==all.n)
    missing.index <- cbind(missing.index,geno[,dat.index])
  }
  colnames(missing.index) <- n.c2
  nmiss <- missing.index + 1
  newSNP <- Impute(nmiss,impute.method ="random")-1
  newZ <- c()
  for(i in 1:nsnp){
    tmp.1 <- unlist(strsplit(SNPna[i],"-"))
    if(length(tmp.1)==1){
      if(substr(tmp.1,1,1)=="a"){
        index1 <- which(gsub("d","",gsub("a","",tmp.1))==n.c2)
        newZ <- cbind(newZ,newSNP[,index1])
      }
      if(substr(tmp.1,1,1)=="d"){
        index1 <- which(gsub("d","",gsub("a","",tmp.1))==n.c2)
        newZ <- cbind(newZ,1-abs(newSNP[,index1]))
      }
    }
    
    if(length(tmp.1)==2){
      if(substr(tmp.1[1],1,1)=="a"){
        index1 <- which(gsub("d","",gsub("a","",tmp.1[1]))==n.c2)
        newZ1 <- newSNP[,index1]
      }
      if(substr(tmp.1[1],1,1)=="d"){
        index1 <- which(gsub("d","",gsub("a","",tmp.1[1]))==n.c2)
        newZ1 <- 1-abs(newSNP[,index1])
      }
      if(substr(tmp.1[2],1,1)=="a"){
        index2 <- which(gsub("d","",gsub("a","",tmp.1[2]))==n.c2)
        newZ2 <- newSNP[,index2]
      }
      if(substr(tmp.1[2],1,1)=="d"){
        index2 <- which(gsub("d","",gsub("a","",tmp.1[2]))==n.c2)
        newZ2 <- 1-abs(newSNP[,index2])
      }
      newZ <- cbind(newZ,newZ1*newZ2)
    }
    
  }
  colnames(newZ) <- SNPna  
  return(newZ)
}


















