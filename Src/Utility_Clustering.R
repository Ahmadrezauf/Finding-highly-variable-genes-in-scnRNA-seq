visualizetSNEggPlot <- function(tsneObj , dims , cellTypes , sampNames , scaleColor , leg){
  res <- 0
  tsneNewCoords <- tsneObj$Y
  if(dims == 2){
    tsneNewCoords <- data.frame(newX = tsneNewCoords[,1] , newY = tsneNewCoords[,2])
    tsneNewCoords$sampNames <- sampNames
    tsneNewCoords$cellType <- cellTypes
    res <- ggplot(tsneNewCoords) + 
      geom_point(aes(x = newX , y = newY , color = cellType)) + 
      scale_color_manual(values = scaleColor)
    if(leg == F){
      res <- res + theme(legend.position = "none")
    }
  }
  
  if(dims == 3){
    require(plotly)
  
    tsneNewCoords <- data.frame(newX = tsneNewCoords[,1] , newY = tsneNewCoords[,2] , newZ = tsneNewCoords[,3])
    tsneNewCoords$sampNames <- sampNames
    tsneNewCoords$cellType <- cellTypes
    
    res <-plot_ly(data = tsneNewCoords , x =~newX , y = ~newY , z = ~newZ , color = ~cellType)
  } 
  res
}

treeStructure <- data.frame()
#resPlots <- list()
vertexLexels <- data.frame(vertex = 0 , level = 0)
fullBinaryTreeStructure <- data.frame()
pieDf <- data.frame()
cellSet <- vector()
nodesID <- 1

resetGlobalVariables <- function(iterations , cellTypes){
  nodesID <<- 1
  treeStructure <<- data.frame(parent = 0 , child = 1)
  fullBinaryTreeStructure <<- data.frame(vertex = 0 , level = 0)
  
  #tmp <- data.frame(x = 1 , y = 1 , name = "Finished")
  #p <- ggplot(tmp) + geom_text(aes(x = x , y = y , label = name , size = x)) + scale_size(range=c(20 , 30))
  # for(i in 1:(2^(iterations + 1) - 1)){
  #   resPlots[[i]] <<- ggplot()
  # }
  pieDf <<- data.frame()
  cellSet <<- unique(cellTypes)
  
  vertexLexels <<- data.frame(vertex = 0 , level = 0)
}

makeVertexLevelTable <- function(iterations , parNodeID , totIters){
  thisNodeID <- nodesID
  nodesID <<- nodesID + 1
  
  
  level <- totIters - iterations + 1
  tmp <- c(thisNodeID , level)
  fullBinaryTreeStructure <<- rbind(fullBinaryTreeStructure , tmp)
  
  if(iterations == 0){
    return(1)
  }
  
  makeVertexLevelTable(iterations = iterations - 1 , parNodeID = thisNodeID , totIters)
  makeVertexLevelTable(iterations = iterations - 1 , parNodeID = thisNodeID , totIters)
}


findHighlyVariableGenes <- function(iData , thresholdForNoSelectedGenes , thisNodeID){
  resultHVG <- HVG_Poisson(iData)
  
  tmp <- resultHVG[[2]]
  
  resultHVG[[2]] %>% 
    group_by(names) %>% 
    mutate(isHVG = ifelse(logAlpha > thresholdForNoSelectedGenes , 1 , 0)) %>% 
    filter(isHVG == 1) -> tmp1
  
  hvgGenes.name <- as.character(tmp1$names)
  idx <- which(rownames(iData) %in% hvgGenes.name)
  
  #write(x = rownames(iData)[idx] , file = paste0("./Results/hvgNames", (nodesID - 1) , ".csv"), append = F)
  
  iData.selected <- iData[idx ,]
  print(paste("Number of selected genes for node" , thisNodeID , "is" , length(idx)))
  
  result <- list(resultHVG[[1]] , resultHVG[[2]] , iData.selected)
}


hvgCluster <- function(iData , iterations , totIters , parNodeID , scaleColor , 
                       thresholdForNoSelectedGenes){
  thisNodeID <- nodesID
  nodesID <<- nodesID + 1
  
  level <- totIters - iterations + 1
  tmp <- c(thisNodeID , level)
  vertexLexels <<- rbind(vertexLexels , tmp)
  
  dcf <- data.frame()
  ## Starting to find HVGs
  tmp <- c(parNodeID , thisNodeID)
  if(parNodeID != 0){
    treeStructure <<- rbind(treeStructure , tmp)  
  }
  tmp <- findHighlyVariableGenes(iData , thresholdForNoSelectedGenes = thresholdForNoSelectedGenes , thisNodeID)

  tmp[[2]]$meansRank <- order(tmp[[2]]$logMeans)
  
  p <- tmp[[1]]
  
  print(p + ggtitle(paste("Logarithm of means , In" , totIters - iterations , "Iterations" , "Having",dim(iData)[2] , "Cells" , 
                          "\n and parent ID is" , parNodeID , "\n And the node ID is", thisNodeID)))
  
  p <- ggplot(tmp[[2]],aes(x = meansRank , y = logAlpha)) + 
    geom_hex() + 
    geom_smooth(method='lm' , color = 'red' , formula = y~x) + 
    xlab("ranks") + ylab("Logarithm of alpha")
  
  print(p + ggtitle(paste("Logarithm of mean Ranks, In" , totIters - iterations , "Iterations" , "Having",dim(iData)[2] , "Cells" , 
                          "\n and parent ID is" , parNodeID , "\n And the node ID is", thisNodeID)))
  
  iData.selected <- tmp[[3]]
  
  # Having duplicates in expression and removing them
  cellTypes <- sapply(strsplit(colnames(iData.selected) , "\\.") , '[' , 1)
  cellCount <- length(cellTypes)
  cellCountTable <- as.data.frame(t(as.matrix(table(cellTypes))))
  zeroCellCounts <- cellSet[!cellSet %in% unique(cellTypes)] 
  cellCountTable[,zeroCellCounts] <- 0
  cellCountTable <- cbind(data.frame(index = thisNodeID , cellCount = cellCount) , cellCountTable)
  
  pieDf <<- rbind(pieDf , cellCountTable)
  
  ## Removing the cells that all genes expression is zero
  tmpColSum <- colSums(iData.selected)
  zeroSamples <- colnames(iData.selected)[tmpColSum == 0]
  iData.selected.Zero <- iData.selected[,(tmpColSum == 0)]
  iData.selected <- iData.selected[,(tmpColSum != 0)]
  zeroCellTypes <- cellTypes[tmpColSum == 0]
  cellTypes <- cellTypes[tmpColSum != 0]
  
  dm2 <- DiffusionMap(data = t(as.matrix(iData.selected)) , n_eigs = 2, density_norm = F , distance = "euclidean")
  
  dcf <- data.frame(DC1 = eigenvectors(dm2)[,1] , DC2 = eigenvectors(dm2)[,2] ,
                    name = colnames(iData.selected) , gr = cellTypes)
  if(length(zeroSamples) != 0){
    dcf <- rbind(dcf , data.frame(DC1 = 0 , DC2 = 0 , name = zeroSamples , gr = zeroCellTypes))
  }
  
  cellTypes <- c(cellTypes , zeroCellTypes)
  iData.selected <- cbind(iData.selected , iData.selected.Zero)
  
  p <- ggplot(dcf) + 
    geom_point(aes(x = DC1 , y = DC2 , color = cellTypes)) + 
    scale_color_manual(values = scaleColor) + 
    theme(legend.position = "none")
  
  print(p)
  
  if(iterations == 0){
    print(paste("Iterations finished for node", thisNodeID))
    return(1)
  }
  
  dcf %>% 
    filter(DC1 > 0) -> tmp1
  
  idx1 <- which(colnames(iData) %in% tmp1$name)
  iData.selected.1 <- as.data.frame(iData[,idx1])
  if(length(idx1) == 1){
    colnames(iData.selected.1) <- colnames(iData)[idx1]
    rownames(iData.selected.1) <- rownames(iData)
  }
  
  idx2 <- which(!(colnames(iData) %in% tmp1$name))
  iData.selected.2 <- as.data.frame(iData[,idx2])
  if(length(idx2) == 1){
    colnames(iData.selected.2) <- colnames(iData)[idx2]
    rownames(iData.selected.2) <- rownames(iData)
  }
  
  flag1 <- ifelse(length(idx1) <= 3 , 1, 0)
  flag2 <- ifelse(length(idx2) <= 3 , 1, 0)
  
  if((flag1 == 1) && (flag2 == 0)){
    
    nodesID <<- nodesID + 1
    tmp <- c(thisNodeID + 1, level + 1)
    vertexLexels <<- rbind(vertexLexels + 1 , tmp)
    
    tmp <- c(thisNodeID , thisNodeID + 1)
    if(parNodeID != 0){
      treeStructure <<- rbind(treeStructure , tmp)  
    }
    
    addToPieDf(iData.selected.1 , thisNodeID + 1)
    
    hvgCluster(iData.selected.2 , iterations - 1 , totIters , parNodeID = thisNodeID , scaleColor = scaleColor , 
               thresholdForNoSelectedGenes)
    
    return(1)
  }
  else if(flag1 == 0 && (flag2 == 1)){
    nodesID <<- nodesID + 1
    tmp <- c(thisNodeID + 1, level + 1)
    vertexLexels <<- rbind(vertexLexels + 1 , tmp)
    
    tmp <- c(thisNodeID , thisNodeID + 1)
    if(parNodeID != 0){
      treeStructure <<- rbind(treeStructure , tmp)  
    }
    
    addToPieDf(iData.selected.2 , thisNodeID + 1)
    
    hvgCluster(iData.selected.1 , iterations - 1 , totIters , parNodeID = thisNodeID , scaleColor = scaleColor , 
               thresholdForNoSelectedGenes)
    return(1)
  }
  else if((flag1 == 1) && (flag2 == 1)){
    return(1)
  }
  
  hvgCluster(iData.selected.1 , iterations - 1 , totIters , parNodeID = thisNodeID , scaleColor = scaleColor , 
             thresholdForNoSelectedGenes)
  hvgCluster(iData.selected.2 , iterations - 1 , totIters , parNodeID = thisNodeID , scaleColor = scaleColor , 
             thresholdForNoSelectedGenes)
  
  print(paste("Iterations done! Node ID" , thisNodeID))
  return(1)
}

addToPieDf <- function(iData , thisNodeID){
  cellTypes <- sapply(strsplit(colnames(iData) , "\\.") , '[' , 1)
  cellCount <- length(cellTypes)
  cellCountTable <- as.data.frame(t(as.matrix(table(cellTypes))))
  zeroCellCounts <- cellSet[!cellSet %in% unique(cellTypes)] 
  cellCountTable[,zeroCellCounts] <- 0
  cellCountTable <- cbind(data.frame(index = thisNodeID , cellCount = cellCount) , cellCountTable)
  
  pieDf <<- rbind(pieDf , cellCountTable)
}
