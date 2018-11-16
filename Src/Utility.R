NormalizeL2Sqr <- function(iData){
  require(dplyr)
  
  sqrData <- sqrt(iData)
  tmp <- iData
  tmp[tmp == 0] <- NA
  meanOfOrgMatrix <- mean(as.matrix(tmp) , na.rm = T)
  meanOfSqrMatrix <- mean(as.matrix(sqrt(tmp)) , na.rm = T)
  remove(tmp)
  sqrData <- as.data.frame(sqrData)
  
  iData %>% 
    colSums(na.rm = T) %>% 
    sqrt() -> sampleSum
  
  if(length(which(sampleSum == 0)) != 0){
    sqrData <- sqrData[,-which(sampleSum == 0)] # removing rows that their sum is zero!   
    sampleSum <- sampleSum[-which(sampleSum == 0)]
    warning("Some samples expression values are zero")
  }
  meanOfSampleSum <- mean(sampleSum , na.rm = T)
  
  # Too slow!
  # sqrData <- apply(sqrData , 2 , function(x) x/sqrt(sum(x^2)))
  if(ncol(sqrData) == 0){
    stop("All of the gene's expression in samples is zero") 
  } else{
    tmp <- lapply(1:ncol(sqrData) , function(i) sqrData[,i] <<- (sqrData[,i]/sampleSum[i]))
    remove(tmp)
  }

  sqrData <- sqrData * meanOfSampleSum
  sqrData
}

NormalizeL2 <- function(iData){
  require(dplyr)
  
  sqrData <- iData
  
  sqrData <- as.data.frame(sqrData)
  
  iData <- iData^2
  iData %>% 
    colSums(na.rm = T) %>% 
    sqrt() -> sampleSum
  
  if(length(which(sampleSum == 0)) != 0){
    sqrData <- sqrData[,-which(sampleSum == 0)] # removing rows that their sum is zero!   
    sampleSum <- sampleSum[-which(sampleSum == 0)]
    warning("Some samples expression values are zero")
  }
  meanOfSampleSum <- mean(sampleSum , na.rm = T)
  
  # Too slow!
  # sqrData <- apply(sqrData , 2 , function(x) x/sqrt(sum(x^2)))
  if(ncol(sqrData) == 0){
    stop("All of the gene's expression in samples is zero") 
  } else{
    tmp <- lapply(1:ncol(sqrData) , function(i) sqrData[,i] <<- (sqrData[,i]/sampleSum[i]))
    remove(tmp)
    sqrData
  }
  sqrData <- sqrData * meanOfSampleSum
  sqrData
}

NormalizeL1 <- function(iData){
  require(dplyr)
  
  sqrData <- iData
  
  sqrData <- as.data.frame(sqrData)
  
  iData %>% 
    colSums(na.rm = T) -> sampleSum
  
  if(length(which(sampleSum == 0)) != 0){
    sqrData <- sqrData[,-which(sampleSum == 0)] # removing rows that their sum is zero!   
    sampleSum <- sampleSum[-which(sampleSum == 0)]
    warning("Some samples expression values are zero")
  }
  meanOfSampleSum <- mean(sampleSum , na.rm = T)
  
  # Too slow!
  # sqrData <- apply(sqrData , 2 , function(x) x/sqrt(sum(x^2)))
  if(ncol(sqrData) == 0){
    stop("All of the gene's expression in samples is zero") 
  }
  else{
    tmp <- lapply(1:ncol(sqrData) , function(i) sqrData[,i] <<- (sqrData[,i]/sampleSum[i]))
    remove(tmp)
    sqrData
  }
  #sqrData <- sqrData * meanOfSampleSum
  sqrData
}

NormalizeL2SqrErcc <- function(iData){
  require(dplyr)
  
  #### WARNING : This is a new version of normalization with l2Square
  
  geneNames <- rownames(iData)
  erccGenes <- geneNames[grep(pattern = "^ERCC" , x = geneNames)]
  notErccRows <- setdiff(x = 1:nrow(iData) , y = which(rownames(iData) %in% erccGenes))
  iDataWithoutErc <- iData[notErccRows,]
  
  iDataWithoutErc <- NormalizeL2Sqr(iDataWithoutErc)
  
  iDataErccGenes <- iData[erccGenes,]
  iDataErccGenes <- NormalizeL2Sqr(iDataErccGenes)
  
  erccSamplesMean <- colMeans(iDataErccGenes , na.rm = T)
  
  #magnifier <- sqrt(nrow(iDataWithoutErc) / nrow(iDataErccGenes))
  #iDataWithoutErc <- iDataWithoutErc * magnifier
  
  iData <- rbind(iDataWithoutErc , iDataErccGenes)
  result <- iData
}

NormalizeL2Ercc <- function(iData){
  require(dplyr)
  
  #### WARNING : This is a new version of normalization with l2Square
  
  geneNames <- rownames(iData)
  erccGenes <- geneNames[grep(pattern = "^ERCC" , x = geneNames)]
  notErccRows <- setdiff(x = 1:nrow(iData) , y = which(rownames(iData) %in% erccGenes))
  iDataWithoutErc <- iData[notErccRows,]
  
  iDataWithoutErc <- NormalizeL2(iDataWithoutErc)
  
  iDataErccGenes <- iData[erccGenes,]
  iDataErccGenes <- NormalizeL2(iDataErccGenes)
  
  erccSamplesMean <- colMeans(iDataErccGenes , na.rm = T)
  
  #magnifier <- sqrt(nrow(iDataWithoutErc) / nrow(iDataErccGenes))
  #iDataWithoutErc <- iDataWithoutErc * magnifier
  
  iData <- rbind(iDataWithoutErc , iDataErccGenes)
  result <- iData
}

flatnessPlotWithErcc <- function(genesInfo , chrTitle){
  genesInfo <- as.data.frame(genesInfo)
  genesInfo$names <- as.character(genesInfo$names)
  erccIndexes <- grep(pattern = "^ERCC" , x = genesInfo$names)
  erccGenes <- genesInfo[erccIndexes , "names"]
  genesInfo %>% 
    mutate(isErcc = ifelse(names %in% erccGenes , 1 , 0)) -> genesInfo
  
  erccGenesInfo <- genesInfo %>% filter(isErcc == 1)
  fit <- lm(formula = logAlpha ~ logMeans , data = erccGenesInfo , na.action = na.omit)
  coefs <- coef(fit)
  
  genesInfoWithERCC <- genesInfo %>% filter(isErcc == 1)
  genesInfo <- genesInfo %>% filter(isErcc == 0)
  p <- ggplot() + 
    geom_point(data = genesInfo , aes(x = logMeans , y = logAlpha) ,  col = "darkslategray3" , alpha = 0.15) +
    geom_point(data = genesInfoWithERCC , aes(x = logMeans , y = logAlpha), , col = "firebrick1") +
    #geom_abline(slope= coefs[2], intercept=coefs[1]) + 
    xlab("Logarithm of means for a gene") + ylab("Logarithm of alpha") + 
    ggtitle(chrTitle)  
}

HVG_Poisson <- function(iData){
  sqrData <- iData
  geneMeans <- rowMeans(sqrData)
  geneVars <- oompaBase::matrixVar(sqrData , geneMeans)
  
  # geneVars <- RowVar(sqrData) # To slow!
  # Also applying a function was tooooo slow!
  geneNames = rownames(sqrData)
  
  remove(sqrData)
  statsDF <- data.frame(means = geneMeans , var = geneVars , names = as.character(geneNames))
  statsDF$logMeans <- log2(statsDF$means)
  statsDF %>% 
    group_by(names) %>% 
    mutate(myAlpha = (var / means)) -> result
  
  result$logAlpha <- log(result$myAlpha)
  
  p <- ggplot(result,aes(x = logMeans , y = logAlpha)) + geom_point(alpha = 0.1) + 
    xlab("Logarithm of gene's average") + ylab("Logarithm of alpha(RV)") + 
    geom_smooth(method = "lm" , formula = y ~ x)

  result <- list(p , result)
  result
}

HVG_NB <- function(iData){
  sqrData <- iData
  geneMeans <- rowMeans(sqrData)
  geneVars <- oompaBase::matrixVar(sqrData , geneMeans)
  
  # geneVars <- RowVar(sqrData) # To slow!
  # Also applying a function was tooooo slow!
  geneNames = rownames(sqrData)
  
  remove(sqrData)
  statsDF <- data.frame(means = geneMeans , var = geneVars , names = as.character(geneNames))
  statsDF$logMeans <- log2(statsDF$means)
  statsDF %>% 
    group_by(names) %>% 
    mutate(myAlpha = ((var - means)/ means^2)) -> result
  
  result$logAlpha <- log(result$myAlpha)
  
  p <- ggplot(result,aes(x = logMeans , y = logAlpha)) + geom_point() + 
    xlab("Logarithm of means for a gene") + ylab("Logarithm of alpha") + 
    geom_smooth(method = "lm" , formula = y ~ x)
  result <- list(p , result)
  result
}

HVG_PoissonOmmitZeros <- function(iData){
  sqrData <- iData
  sqrData[sqrData == 0] <- NA
  geneMeans <- rowMeans(sqrData , na.rm = T)
  geneVars <- oompaBase::matrixVar(sqrData , geneMeans , na.rm = T)
  
  # geneVars <- RowVar(sqrData) # To slow!
  # Also applying a function was tooooo slow!
  geneNames = rownames(sqrData)
  
  remove(sqrData)
  statsDF <- data.frame(means = geneMeans , var = geneVars , names = geneNames)
  statsDF$logMeans <- log2(statsDF$means)
  statsDF %>% 
    group_by(names) %>% 
    mutate(myAlpha = (var / means^2)) -> result
  
  result$logAlpha <- log(result$myAlpha)
  
  p <- ggplot(result,aes(x = logMeans , y = logAlpha)) + geom_point(alpha = 0.1) + 
    xlab("Logarithm of gene's average") + ylab("Logarithm of CV2") + 
    geom_smooth(method = "lm" , formula = y ~ x)
  
  result <- list(p , result)
  result
}

HVG_NBOmmitZeros <- function(iData){
  
  sqrData <- iData
  sqrData[sqrData == 0] <- NA
  geneMeans <- rowMeans(sqrData , na.rm = T)
  geneVars <- oompaBase::matrixVar(sqrData , geneMeans , na.rm = T)
  
  # geneVars <- RowVar(sqrData) # To slow!
  # Also applying a function was tooooo slow!
  geneNames = rownames(sqrData)
  
  remove(sqrData)
  statsDF <- data.frame(means = geneMeans , var = geneVars , names = geneNames)
  statsDF$logMeans <- log2(statsDF$means)
  statsDF %>% 
    group_by(names) %>% 
    mutate(myAlpha = ((var - means) / means^2)) -> result
  
  result$logAlpha <- log(result$myAlpha)
  
  p <- ggplot(result,aes(x = logMeans , y = logAlpha)) + geom_point(alpha = 0.2) + 
    xlab("Logarithm of gene's average") + ylab("Logarithm of alpha") + 
    geom_smooth(method = "lm" , formula = y ~ x)
  
  result <- list(p , result)
  result
}

HVG_residual <- function(iData){
  require(ggplot2)
  require(dplyr)
  require(ggplot2)
  require(TailRank)
  
  sqrData <- iData
  
  geneMeans <- rowMeans(sqrData , na.rm = T)
  geneVars <- oompaBase::matrixVar(sqrData , geneMeans , na.rm = T)
  
  geneNames = rownames(sqrData)
  
  remove(sqrData)
  statsDF <- data.frame(means = geneMeans , var = geneVars , names = geneNames)
  
  statsDF %>% 
    filter(means != 0) %>% 
    filter(var != 0) -> statsDF
  statsDF$logMeans <- log2(statsDF$means)
  statsDF$logVars <- log2(statsDF$var)
  result <- statsDF
  
  ## testing for log means and log vars
  fit <- lm(logVars ~ logMeans, data = result , na.action = na.omit)
  coefs <- coef(fit)
  summary(fit)
  
  resids <- residuals(fit)
  result$resids <- resids
  result %>% group_by(names) %>% 
    mutate(isAbove = ifelse(logVars > (coefs[1] + coefs[2] * logMeans) , 1 , 0)) -> result

  p <- ggplot(result,aes(x = logMeans , y = logVars)) + 
    geom_point(alpha = 0.05) + 
    geom_smooth(formula = y~x , color = "red" , method = "lm") + 
    xlab("Logarithm of gene's expression average") + 
    ylab("Logarithm of gene's expression estimated variance")
  
  p <- ggplot(result) + 
    geom_point(aes(x = logMeans , y = resids) , alpha = 0.05) +
    xlab("Logarithm of gene's average") + 
    ylab("Residual from regression line") + 
    ggtitle("Flatness of residuals for regression of logmeans and logvars")

  # ggplot(result) + 
  #   geom_hex(aes(x = logMeans , y = resids))
  res <- list(p, result)
  res
}


HVG_residualTwoTerms <- function(iData , signal){
  require(ggplot2)
  require(dplyr)
  require(ggplot2)
  require(TailRank)
  
  sqrData <- iData
  
  geneMeans <- rowMeans(sqrData , na.rm = T)
  geneVars <- oompaBase::matrixVar(sqrData , geneMeans , na.rm = T)
  
  geneNames = rownames(sqrData)
  
  remove(sqrData)
  statsDF <- data.frame(means = geneMeans , var = geneVars , names = geneNames)
  
  statsDF %>% 
    filter(means != 0) %>% 
    filter(var != 0) -> statsDF
  statsDF$logMeans <- log2(statsDF$means)
  statsDF$logVars <- log2(statsDF$var)
  result <- statsDF
  
  ## testing for log means and log vars
  fit <- lm(logVars ~ (logMeans + I(logMeans^2))  , data = result , na.action = na.omit)
  coefs <- coef(fit)
  summary(fit)
  
  if(signal == 1){
    p <- ggplot(result,aes(x = logMeans , y = logVars)) + 
      geom_point() + 
      geom_smooth(method = 'lm' , color = 'red' , formula = y~(x + I(x^2)) , na.rm = T) + 
      xlab("Means for a gene") + ylab("Variance")
    print(p)
  }
  
  resids <- residuals(fit)
  result$resids <- resids
  result %>% group_by(names) %>% 
    mutate(isAbove = ifelse(logVars > (coefs[1] + coefs[2] * logMeans + coefs[3] * logMeans^2) , 1 , 0)) -> result
  
  p <- ggplot(result) + 
    geom_point(aes(x = logMeans , y = resids)) + 
    ggtitle("Flatness of residuals for regression of log of means and log of variances")
  
  resultTmp <- result
  colnames(result)[which(colnames(result) == "resids")] <- "logAlpha"
  
  #f <- flatnessPlotWithErcc(result , "Flatness")
  # ggplot(result) + 
  #   geom_hex(aes(x = logMeans , y = resids))
  res <- list(p, result)
  res
}

